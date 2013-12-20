#include <cuda.h>
#include <stdio.h>
#include <cufft.h>
#include <string.h>
#include "utils_cuda.h"
//#include <cutil_inline.h>
#include "data.h"
#include "timer.h"


__global__ void Fir_SpB_wR_Sep_Shared(float *d_data, float *d_coeff, int nTaps, int nChannels, float *d_spectra) {
	int c = 0;
	int ypos = blockDim.x*blockIdx.y;
	short int ccounter=257,tcounter=257;
	short int column,oldesttap;
	float ftemp=0;
	int bl=0;
	if(blockIdx.x!=0) bl=blockIdx.x*nChannels*125+threadIdx.x+ypos-7*nChannels;
	__shared__ float l_buffer[128*8];
	//__shared__ float l_coeff[128*8];
	l_buffer[threadIdx.x]=d_data[bl];
	l_buffer[128+threadIdx.x]=d_data[bl+nChannels];
	l_buffer[256+threadIdx.x]=d_data[bl+nChannels*2];
	l_buffer[384+threadIdx.x]=d_data[bl+nChannels*3];
	l_buffer[512+threadIdx.x]=d_data[bl+nChannels*4];
	l_buffer[640+threadIdx.x]=d_data[bl+nChannels*5];
	l_buffer[768+threadIdx.x]=d_data[bl+nChannels*6];
	oldesttap=7;
	
	//load buffer
	__syncthreads();
	//if (threadIdx.x==0 && blockIdx.x==0) printf("test");
	
	for(int bl=(blockIdx.x*nChannels*125+threadIdx.x+ypos);bl<(blockIdx.x+125)*nChannels;bl+=nChannels){
		//update buffer
		//l_buffer[oldesttap*128+threadIdx.x]=d_data[bl+threadIdx.x];
		l_buffer[threadIdx.x]=d_data[bl];
		__syncthreads();
		oldesttap++;
		tcounter=tcounter<<1;
		//oldesttap=oldesttap*((tcounter|1)^1);
		oldesttap=oldesttap%8;
		ftemp=0.0f;
		column=oldesttap;
		for(int t=0;t<nTaps;t++){
			c=column*128 + threadIdx.x;
			//ftemp  += d_coeff[column*nChannels + ypos + threadIdx.x]*l_buffer[c];
			ftemp  += d_coeff[column*nChannels + threadIdx.x]*l_buffer[c];
			column++;
			ccounter=ccounter<<1;
			//column=column*((ccounter|1)^1);
			column=column%8;
		}
		
		c=bl + threadIdx.x + ypos;
		d_spectra[c]=ftemp;
	}
}

__global__ void Fir(float *d_signal_real, float *d_signal_img, const float* coeff, const int nTaps, const int nChannels, float2 *spectra)
{
	int tx = threadIdx.x + blockDim.x*blockIdx.y;
	int index = nChannels*blockIdx.x + tx;
	int i, i_coeff, i_data;
	float local_spectra_x = 0.0f;
	float local_spectra_y = 0.0f;
	
	for(int t=0;t<nTaps;t++){
	  i = t*nChannels;
	  i_coeff = i + tx;
	  i_data = i + index;
	  local_spectra_x += coeff[i_coeff]*d_signal_real[i_data];
	  local_spectra_y += coeff[i_coeff]*d_signal_img[i_data];
	}
		
	spectra[index].x=local_spectra_x;
	spectra[index].y=local_spectra_y;
	//return;
}

__global__ void Fir_01(float *d_signal_real, float *d_signal_img, const float* coeff, const int nTaps, const int nChannels, float2 *spectra, unsigned int datasize)
{
	int tx = threadIdx.x + blockDim.x*blockIdx.y;
	int i, i_coeff, i_data;
	float local_spectra_x = 0.0f;
	float local_spectra_y = 0.0f;
 
	for (int index = nChannels*blockIdx.x + tx; index < datasize - nTaps*nChannels; index+=gridDim.x*nChannels){
	for(int t=0;t<nTaps;t++){
	  i = t*nChannels;
	  i_coeff = i + tx;
	  i_data = i + index;
	  local_spectra_x += coeff[i_coeff]*d_signal_real[i_data];
	  local_spectra_y += coeff[i_coeff]*d_signal_img[i_data];
	}
		
	spectra[index].x=local_spectra_x;
	spectra[index].y=local_spectra_y;
	}
	return;
}

void gpu_code(  float *real,
	        float *img, 
	        float2 *spectra, 
		float *coeff,
		const int nChannels,
		unsigned int nBlocks, 
		unsigned int filesize,
		int blocks_y){
//------------ initialize card -----------

  int devCount, device;
  
  checkCudaErrors(cudaGetDeviceCount(&devCount));
  printf("\n\t\t-------------- GPU part -----------------");
  printf("\nThere are %d devices.", devCount);

  for (int i = 0; i < devCount - 1; i++){
	cudaDeviceProp devProp;
	checkCudaErrors(cudaGetDeviceProperties(&devProp,i));	
	printf("\n\t Using device:\t\t\t%s\n", devProp.name);
	device = i;
	checkCudaErrors(cudaSetDevice(device));
  }

//------------ memory setup -------------------------------------
	float2 *d_spectra;
	float  *d_real, *d_img, *d_coeff, *d_spectra2;
	GpuTimer timer;

	float run_time = -1.1f;

	//malloc
	printf("\nDevice memory allocation...\t\t");
	timer.Start();
	checkCudaErrors(cudaMalloc((void **) &d_spectra, sizeof(float2)*filesize));
	checkCudaErrors(cudaMalloc((void **) &d_spectra2, 2*sizeof(float)*filesize));
	checkCudaErrors(cudaMalloc((void **) &d_coeff,   sizeof(float)*nChannels*nTaps));
	checkCudaErrors(cudaMalloc((void **) &d_real,    sizeof(float)*filesize));
	checkCudaErrors(cudaMalloc((void **) &d_img,     sizeof(float)*filesize));
	timer.Stop();
	printf("done in %g ms.", timer.Elapsed());

	// set to 0.0
	printf("\nDevice memset...\t\t\t");
	timer.Start();
	checkCudaErrors(cudaMemset(d_spectra, 0.0, sizeof(float2)*filesize));
	checkCudaErrors(cudaMemset(d_spectra2, 0.0, 2*sizeof(float)*filesize));
	timer.Stop();
	printf("done in %g ms.", timer.Elapsed());

	// copy data to device
	printf("\nCopy data from host to device...\t");
	timer.Start();
	checkCudaErrors(cudaMemcpy(d_coeff, coeff, nChannels*nTaps*sizeof(float), cudaMemcpyHostToDevice));
	//checkCudaErrors(cudaMemcpyToSymbol(C_coeff,  coeff,   sizeof(float)*nChannels*nTaps));
	checkCudaErrors(cudaMemcpy(d_real, real, filesize*sizeof(float), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_img,  img,  filesize*sizeof(float), cudaMemcpyHostToDevice));
	timer.Stop();
	printf("done in %g ms.", timer.Elapsed());
//---------------------------------------------------------

//--------------- Fir ----------------------------
	//dim3 gridSize(nBlocks - nTaps + 1, blocks_y, 1);
	dim3 gridSize(480, 8, 1);
	//dim3 gridSize(512, blocks_y, 1);
	//dim3 blockSize(nChannels/gridSize.y, 1, 1); 
	dim3 blockSize(128, 1, 1); 
	
	timer.Start();
	Fir<<<gridSize, blockSize>>>(d_real, d_img, d_coeff, nTaps, nChannels, d_spectra);
	timer.Stop();
	run_time=timer.Elapsed();
	cudaDeviceSynchronize(); checkCudaErrors(cudaGetLastError());
	printf("\n\t\t------------ Kernel run -----------------");
	printf("\nFir kernel \n");
	printf("\n\n blocks \t time \t\t threads \t bandwidth \t flops");
	printf("\n%d \t\t %lf \t %i \t\t %g \t %g\n",nBlocks/12,run_time,nChannels/blocks_y, 53248.0*(nBlocks-nTaps+1)*1000/run_time, 16384.0*(nBlocks-nTaps+1)*1000.0/run_time);
	Fir_SpB_wR_Sep_Shared<<<gridSize, blockSize>>>(d_real, d_coeff, nTaps, nChannels, d_spectra2);

//--------------- cuFFT ----------------------------
/*
	//Create fft Plan
	cufftHandle plan;
	cufftPlan1d(&plan, nChannels, CUFFT_C2C, nBlocks);

	//execute plan and copy back to host
	printf("\n\ncuFFT plan...\t\t");
	timer.Start();
	cufftExecC2C(plan, (cufftComplex *)d_spectra, (cufftComplex *)d_spectra, CUFFT_FORWARD);
	timer.Stop();
	printf("done in %g ms.\n\n", timer.Elapsed());

	//Destroy the cuFFT plan
	cufftDestroy(plan);
*/
//--------------- copy data back ----------------------------

	//checkCudaErrors(cudaMemcpy(spectra,d_spectra,filesize*sizeof(float2), cudaMemcpyDeviceToHost));	
	checkCudaErrors(cudaMemcpy(spectra,d_spectra2,2*filesize*sizeof(float), cudaMemcpyDeviceToHost));	

//--------------- clean-up process ----------------------------
	checkCudaErrors(cudaFree(d_spectra));
	checkCudaErrors(cudaFree(d_spectra2));
	checkCudaErrors(cudaFree(d_real));
	checkCudaErrors(cudaFree(d_img));
	checkCudaErrors(cudaFree(d_coeff));

}
