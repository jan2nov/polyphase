#include <cuda.h>
#include <stdio.h>
#include <cufft.h>
#include <string.h>
#include "utils_cuda.h"
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
  printf("There are %d devices.", devCount);

  for (int i = 0; i < devCount; i++){
    cudaDeviceProp devProp;
    checkCudaErrors(cudaGetDeviceProperties(&devProp, i));
    if (!strncmp("Tesla", devProp.name, 5)) device = i;
  }
  printf("\n\t Using device:\t\t\t%d", device);

  checkCudaErrors(cudaSetDevice(device));
//------------ memory setup -------------------------------------
	float2 *d_spectra;
	float  *d_real, *d_img, *d_coeff;

	float run_time = -1.1f;

	//malloc
	checkCudaErrors(cudaMalloc((void **) &d_spectra, sizeof(float2)*filesize));
	checkCudaErrors(cudaMalloc((void **) &d_coeff,   sizeof(float)*nChannels*nTaps));
	checkCudaErrors(cudaMalloc((void **) &d_real,    sizeof(float)*filesize));
	checkCudaErrors(cudaMalloc((void **) &d_img,     sizeof(float)*filesize));

	// set to 0.0
	checkCudaErrors(cudaMemset(d_spectra, 0.0, sizeof(float2)*filesize));

	// copy data to device
	checkCudaErrors(cudaMemcpy(d_coeff, coeff, nChannels*nTaps*sizeof(float), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_real, real, filesize*sizeof(float), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_img,  img,  filesize*sizeof(float), cudaMemcpyHostToDevice));

//---------------------------------------------------------

//--------------- Fir ----------------------------
	dim3 gridSize(nBlocks - nTaps + 1, blocks_y, 1);
	//dim3 gridSize(512, blocks_y, 1);
	dim3 blockSize(nChannels/gridSize.y, 1, 1); 
	
	GpuTimer timer;

	timer.Start();
	Fir<<<gridSize, blockSize>>>(d_real, d_img, d_coeff, nTaps, nChannels, d_spectra);
	timer.Stop();
	run_time=timer.Elapsed();
	cudaDeviceSynchronize(); checkCudaErrors(cudaGetLastError());

	printf("%d %lf %i %lf %lf\n",nBlocks/12,run_time,nChannels/blocks_y, 53248.0*(nBlocks-nTaps+1)*1000/run_time, 16384.0*(nBlocks-nTaps+1)*1000.0/run_time);

//--------------- cuFFT ----------------------------

	//Create fft Plan
	cufftHandle plan;
	cufftPlan1d(&plan, nChannels, CUFFT_C2C, nBlocks);

	//execute plan and copy back to host
	cufftExecC2C(plan, (cufftComplex *)d_spectra, (cufftComplex *)d_spectra, CUFFT_FORWARD);
	

	//Destroy the cuFFT plan
	cufftDestroy(plan);

//--------------- copy data back ----------------------------

	checkCudaErrors(cudaMemcpy(spectra,d_spectra,filesize*sizeof(float2), cudaMemcpyDeviceToHost));	

//--------------- clean-up process ----------------------------
	checkCudaErrors(cudaFree(d_spectra));
	checkCudaErrors(cudaFree(d_real));
	checkCudaErrors(cudaFree(d_img));
	checkCudaErrors(cudaFree(d_coeff));

}
