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
  GpuTimer timer, time_memory, time_kernels;
  
  checkCudaErrors(cudaGetDeviceCount(&devCount));
  printf("\n\t\t-------------- GPU part -----------------");
  printf("\nThere are %d devices.", devCount);

	//get number of GPU available
  for (int i = 0; i < devCount; i++){
	cudaDeviceProp devProp;
	checkCudaErrors(cudaGetDeviceProperties(&devProp,i));	
	printf("\n\t Using device:\t\t\t%s\n", devProp.name);
	device = i;
	// set some preferable card
	checkCudaErrors(cudaSetDevice(device));
  }


//------------- stream setup ------------------------------------
timer.Start();
	cudaStream_t stream0, stream1;
 
	printf("\nStream creating...\t\t\t");
	timer.Start();
		checkCudaErrors(cudaStreamCreate(&stream0));
		checkCudaErrors(cudaStreamCreate(&stream1));
	timer.Stop();
	printf("done in %g ms.", timer.Elapsed());
	

//---------------------------------------------------------------

//------------ memory setup -------------------------------------
	float *d_coeff;
	int seg_blocks = 5000; // each segment compute # of spectra
	int run_blocks = nBlocks - nTaps + 1; // needed blocks for run on whole host data
	int SegSize = (seg_blocks + nTaps - 1)*nChannels; //size of each segment in the buffer
	int seg_offset = seg_blocks*nChannels;
	
	//stream 0
	float  *d_real_0, *d_img_0; 
	float2 *d_spectra_0;
	
	//stream 1	
	float  *d_real_1, *d_img_1; 
	float2 *d_spectra_1;

	float run_time = 0.0f;
	float mem_time = 0.0f;
	float ker_time = 0.0f;

	// grid and block size
	//int grid_0 = (int)(seg_blocks/2);
	dim3 gridSize0( seg_blocks, blocks_y, 1);
	dim3 blockSize0(nChannels/gridSize0.y, 1, 1); 
	//dim3 gridSize1( run_blocks - grid_0, blocks_y, 1);
	//dim3 blockSize1(nChannels/gridSize1.y, 1, 1); 
		
	checkCudaErrors(cudaMalloc((void **) &d_spectra_0, sizeof(float2)*SegSize));
	checkCudaErrors(cudaMalloc((void **) &d_spectra_1, sizeof(float2)*SegSize));
	checkCudaErrors(cudaMalloc((void **) &d_coeff,   sizeof(float)*nChannels*nTaps));
	checkCudaErrors(cudaMalloc((void **) &d_real_0,    sizeof(float)*SegSize));
	checkCudaErrors(cudaMalloc((void **) &d_img_0,     sizeof(float)*SegSize));
	checkCudaErrors(cudaMalloc((void **) &d_real_1,    sizeof(float)*SegSize));
	checkCudaErrors(cudaMalloc((void **) &d_img_1,     sizeof(float)*SegSize));
	
	//coefficients copy
	checkCudaErrors(cudaMemcpy(d_coeff, coeff, nChannels*nTaps*sizeof(float), cudaMemcpyHostToDevice));

	// set to 0.0

	checkCudaErrors(cudaMemset(d_spectra_0, 0.0, sizeof(float2)*SegSize));
	checkCudaErrors(cudaMemset(d_spectra_1, 0.0, sizeof(float2)*SegSize));

	// copy data to device
for (int i = 0; i < run_blocks; i+=seg_blocks*2){
	
	time_memory.Start();
		checkCudaErrors(cudaMemcpyAsync(d_real_0, real + i*nChannels, sizeof(float)*SegSize, cudaMemcpyHostToDevice, stream0));
		checkCudaErrors(cudaMemcpyAsync(d_img_0, img + i*nChannels, sizeof(float)*SegSize, cudaMemcpyHostToDevice, stream0));
		checkCudaErrors(cudaMemcpyAsync(d_real_1, real + seg_offset + i*nChannels, sizeof(float)*SegSize, cudaMemcpyHostToDevice, stream1));
		checkCudaErrors(cudaMemcpyAsync(d_img_1, img + seg_offset + i*nChannels, sizeof(float)*SegSize, cudaMemcpyHostToDevice, stream1));
	time_memory.Stop();
	mem_time+=time_memory.Elapsed();

	time_kernels.Start();
		Fir<<<gridSize0, blockSize0, 0, stream0>>>(d_real_0, d_img_0, d_coeff, nTaps, nChannels, d_spectra_0);
		Fir<<<gridSize0, blockSize0, 0, stream1>>>(d_real_1, d_img_1, d_coeff, nTaps, nChannels, d_spectra_1);
	time_kernels.Stop();
	ker_time+=time_kernels.Elapsed();

	time_memory.Start();
		checkCudaErrors(cudaMemcpyAsync(spectra + i*nChannels, d_spectra_0, sizeof(float2)*SegSize, cudaMemcpyDeviceToHost, stream0));
		checkCudaErrors(cudaMemcpyAsync(spectra + i*nChannels + seg_offset, d_spectra_1, sizeof(float2)*SegSize, cudaMemcpyDeviceToHost, stream1));
	time_memory.Stop();
	mem_time+=time_memory.Elapsed();
	//printf("\nfilesize: %i Run %i\t\n", filesize, i);
}
	timer.Stop();
	run_time=timer.Elapsed();
	printf("\nDone in %g ms.\nMemory in %g ms.\nKernels in %g ms.\n", run_time, mem_time, ker_time);
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

//--------------- clean-up process ----------------------------
	checkCudaErrors(cudaFree(d_spectra_0));
	checkCudaErrors(cudaFree(d_real_0));
	checkCudaErrors(cudaFree(d_img_0));
	checkCudaErrors(cudaFree(d_spectra_1));
	checkCudaErrors(cudaFree(d_real_1));
	checkCudaErrors(cudaFree(d_img_1));
	checkCudaErrors(cudaFree(d_coeff));
	
	checkCudaErrors(cudaStreamDestroy(stream0));
	checkCudaErrors(cudaStreamDestroy(stream1));

}
