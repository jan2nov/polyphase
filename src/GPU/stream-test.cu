#include <cuda.h>
#include <stdio.h>
#include <cufft.h>
#include <string.h>
#include "utils_cuda.h"
#include "utils_file.h"
//#include <cutil_inline.h>
#include "data.h"
#include "timer.h"

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

__global__ void Fir_SpB(float2* d_data, float* d_coeff, int nTaps, int nChannels, int yshift, float2* d_spectra) {
	int t = 0;
	int bl= blockIdx.x*nChannels;
	int ypos = blockDim.x*blockIdx.y + yshift;
	float2 ftemp1;
	ftemp1.x=0.0f;ftemp1.y=0.0f;

	for(t=ypos + threadIdx.x;t<nTaps*nChannels;t+=nChannels){
		ftemp1.x  += d_coeff[t]*d_data[bl+t].x;
		ftemp1.y  += d_coeff[t]*d_data[bl+t].y;
	}

	t=bl + ypos + threadIdx.x;
	d_spectra[t]=ftemp1;
	return;
}


void gpu_code(  float2 *data, 
				float2 *spectra, 
				float *coeff,
				const int nChannels,
				unsigned int nBlocks, 
				unsigned int filesize,
				int blocks_y, 
				int nTaps, 
				int seg_blocks){
					
bool WRITE=true;
//------------ initialize card -----------

  int devCount, device;
  cudaDeviceProp devProp;
    
  checkCudaErrors(cudaGetDeviceCount(&devCount));
  printf("\n\t\t-------------- GPU part -----------------");
  printf("\nThere are %d devices.", devCount);

	//get number of GPU available
  for (int i = 0; i < devCount; i++){
	checkCudaErrors(cudaGetDeviceProperties(&devProp,i));	
	printf("\n\t Using device:\t\t\t%s\n", devProp.name);
	printf("\t Concurrent kernels:\t\t\t%i\n", devProp.concurrentKernels);
	printf("\t Async Engine Count:\t\t\t%i\n", devProp.asyncEngineCount);
	device = 0;
	// set some preferable card
	checkCudaErrors(cudaSetDevice(device));
  }

	GpuTimer timer, time_memory, time_kernels;

//------------- stream setup ------------------------------------
	cudaStream_t stream0, stream1, stream2, stream3;
 
	printf("\nStream creating...\t\t\t");
	timer.Start();
			checkCudaErrors(cudaStreamCreate(&stream0));
			checkCudaErrors(cudaStreamCreate(&stream1));
			checkCudaErrors(cudaStreamCreate(&stream2));
			checkCudaErrors(cudaStreamCreate(&stream3));
	timer.Stop();
	printf("done in %g ms.\n", timer.Elapsed());
	

//---------------------------------------------------------------

//------------ memory setup -------------------------------------
	float *d_coeff;
	//int seg_blocks = 10000; // each segment compute # of spectra
	int run_blocks = nBlocks - nTaps + 1; // needed blocks for run on whole host data
	int SegSize = (seg_blocks + nTaps - 1)*nChannels; //size of each segment in the buffer
	int seg_offset = seg_blocks*nChannels;
	printf("Number of spectra per block: \t%i\n", seg_blocks);
	printf("Size of segment in bytes: \t%lu\n", SegSize*sizeof(float));
	printf("Run_blocks: \t\t\t%i\n", run_blocks);
	printf("Offset: \t\t\t%i\n", seg_offset);
	printf("-----------------------------------------\n");
	
	
	//stream 0..4
	float2 *d_spectra_0, *d_data_0;
	float2 *d_spectra_1, *d_data_1;
	float2 *d_spectra_2, *d_data_2;
	float2 *d_spectra_3, *d_data_3;

	float fir_time = 0.0f;
	float fft_time = 0.0f;
	float mem_time_in = 0.0f;
	float mem_time_out = 0.0f;
	

	// grid and block size
	//int grid_0 = (int)(seg_blocks/2);
	dim3 gridSize0( seg_blocks, blocks_y, 1);
	dim3 blockSize0(nChannels/gridSize0.y, 1, 1); 
	//dim3 gridSize1( run_blocks - grid_0, blocks_y, 1);
	//dim3 blockSize1(nChannels/gridSize1.y, 1, 1); 
		
	checkCudaErrors(cudaMalloc((void **) &d_spectra_0, sizeof(float2)*SegSize));
	checkCudaErrors(cudaMalloc((void **) &d_spectra_1, sizeof(float2)*SegSize));
	checkCudaErrors(cudaMalloc((void **) &d_spectra_2, sizeof(float2)*SegSize));
	checkCudaErrors(cudaMalloc((void **) &d_spectra_3, sizeof(float2)*SegSize));
	checkCudaErrors(cudaMalloc((void **) &d_coeff,   sizeof(float)*nChannels*nTaps));
	checkCudaErrors(cudaMalloc((void **) &d_data_0,    sizeof(float2)*SegSize));	
	checkCudaErrors(cudaMalloc((void **) &d_data_1,    sizeof(float2)*SegSize));
	checkCudaErrors(cudaMalloc((void **) &d_data_2,    sizeof(float2)*SegSize));	
	checkCudaErrors(cudaMalloc((void **) &d_data_3,    sizeof(float2)*SegSize));
	
	printf("\n\t\td_spectra using filesize: \t%g MB.", 4*sizeof(float2)*SegSize/1024.0/1024);
	printf("\n\t\td_coeff using filesize: \t%g MB.", sizeof(float)*SegSize/1024.0/1024);
	printf("\n\t\td_data using filesize: \t\t%g MB.", 4*sizeof(float2)*SegSize/1024.0/1024);
	printf("\n\t\t----------------------\t\t-----------");
	printf("\n\t\tTotal: \t\t\t\t%g MB.\n\n",sizeof(float)*(17.0*SegSize)/1024/1024);
	
	//coefficients copy
	checkCudaErrors(cudaMemcpy(d_coeff, coeff, nChannels*nTaps*sizeof(float), cudaMemcpyHostToDevice));

	// set to 0.0
	//checkCudaErrors(cudaMemset(d_spectra_0, 0.0, sizeof(float2)*SegSize));
	//checkCudaErrors(cudaMemset(d_spectra_1, 0.0, sizeof(float2)*SegSize));

	//Create fft Plan
	cufftHandle plan0;
	cufftHandle plan1;
	cufftHandle plan2;
	cufftHandle plan3;
	cufftPlan1d(&plan0, nChannels, CUFFT_C2C, seg_blocks);
	cufftPlan1d(&plan1, nChannels, CUFFT_C2C, seg_blocks);
	cufftPlan1d(&plan2, nChannels, CUFFT_C2C, seg_blocks);
	cufftPlan1d(&plan3, nChannels, CUFFT_C2C, seg_blocks);
	cufftSetStream(plan0,stream0);
	cufftSetStream(plan1,stream1);
	cufftSetStream(plan2,stream2);
	cufftSetStream(plan3,stream3);
	
	timer.Start();
for (int i = 0; i < run_blocks; i+=seg_blocks*4){

		checkCudaErrors(cudaMemcpyAsync(d_data_0, data + i*nChannels, sizeof(float2)*SegSize, cudaMemcpyHostToDevice, stream0));
		checkCudaErrors(cudaMemcpyAsync(d_data_1, data + seg_offset + i*nChannels, sizeof(float2)*SegSize, cudaMemcpyHostToDevice, stream1));
		checkCudaErrors(cudaMemcpyAsync(d_data_2, data + 2*seg_offset + i*nChannels, sizeof(float2)*SegSize, cudaMemcpyHostToDevice, stream2));
		checkCudaErrors(cudaMemcpyAsync(d_data_3, data + 3*seg_offset + i*nChannels, sizeof(float2)*SegSize, cudaMemcpyHostToDevice, stream3));

		Fir_SpB<<<gridSize0, blockSize0, 0, stream0>>>(d_data_0, d_coeff, nTaps, nChannels, 0, d_spectra_0);
		cufftExecC2C(plan0, (cufftComplex *)d_spectra_0, (cufftComplex *)d_spectra_0, CUFFT_FORWARD);
		checkCudaErrors(cudaMemcpyAsync(spectra + i*nChannels, d_spectra_0, sizeof(float2)*SegSize, cudaMemcpyDeviceToHost, stream0));

		Fir_SpB<<<gridSize0, blockSize0, 0, stream1>>>(d_data_1, d_coeff, nTaps, nChannels, 0, d_spectra_1);
		cufftExecC2C(plan1, (cufftComplex *)d_spectra_1, (cufftComplex *)d_spectra_1, CUFFT_FORWARD);
		checkCudaErrors(cudaMemcpyAsync(spectra + i*nChannels + seg_offset, d_spectra_1, sizeof(float2)*SegSize, cudaMemcpyDeviceToHost, stream1));
		
		Fir_SpB<<<gridSize0, blockSize0, 0, stream2>>>(d_data_2, d_coeff, nTaps, nChannels, 0, d_spectra_2);
		cufftExecC2C(plan2, (cufftComplex *)d_spectra_2, (cufftComplex *)d_spectra_2, CUFFT_FORWARD);
		checkCudaErrors(cudaMemcpyAsync(spectra + i*nChannels + 2*seg_offset, d_spectra_2, sizeof(float2)*SegSize, cudaMemcpyDeviceToHost, stream2));
		
		Fir_SpB<<<gridSize0, blockSize0, 0, stream3>>>(d_data_3, d_coeff, nTaps, nChannels, 0, d_spectra_3);
		cufftExecC2C(plan3, (cufftComplex *)d_spectra_3, (cufftComplex *)d_spectra_3, CUFFT_FORWARD);
		checkCudaErrors(cudaMemcpyAsync(spectra + i*nChannels + 3*seg_offset, d_spectra_3, sizeof(float2)*SegSize, cudaMemcpyDeviceToHost, stream3));
}

	timer.Stop();
	fir_time=timer.Elapsed();
	printf("\nDone in %g ms.\n", fir_time);

//---------------- write to file process ----------------------

	char str[200];
	sprintf(str,"GPU-stream-%s.dat",devProp.name);
	
		printf("\n Write results into file...\t");
		if (WRITE) save_time(str, nBlocks-nTaps+1, fir_time, fft_time, mem_time_in, mem_time_out, nChannels, nTaps);
		printf("\t done.\n-------------------------------------\n");


//--------------- clean-up process ----------------------------
	checkCudaErrors(cudaFree(d_spectra_0));
	checkCudaErrors(cudaFree(d_spectra_1));
	checkCudaErrors(cudaFree(d_spectra_2));
	checkCudaErrors(cudaFree(d_spectra_3));
	checkCudaErrors(cudaFree(d_data_0));
	checkCudaErrors(cudaFree(d_data_1));
	checkCudaErrors(cudaFree(d_data_2));
	checkCudaErrors(cudaFree(d_data_3));
	checkCudaErrors(cudaFree(d_coeff));

	cufftDestroy(plan0);
	cufftDestroy(plan1);
	cufftDestroy(plan2);
	cufftDestroy(plan3);

	checkCudaErrors(cudaStreamDestroy(stream0));
	checkCudaErrors(cudaStreamDestroy(stream1));
	checkCudaErrors(cudaStreamDestroy(stream2));
	checkCudaErrors(cudaStreamDestroy(stream3));

}
