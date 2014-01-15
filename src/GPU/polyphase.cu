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

void gpu_code(  float *real,
				float *img, 
				float2 *spectra, 
				float *coeff,
				const int nChannels,
				unsigned int nBlocks, 
				unsigned int filesize,
				int blocks_y){

 bool WRITE=true;					
//------------ initialize card -----------

  int devCount, device;
  int maxgrid_x = 0;
  
  checkCudaErrors(cudaGetDeviceCount(&devCount));
  printf("\n\t\t-------------- GPU part -----------------");
  printf("\nThere are %d devices.", devCount);

	for (int i = 0; i < devCount-3; i++){
		cudaDeviceProp devProp;
		device = 0;
		checkCudaErrors(cudaGetDeviceProperties(&devProp,device));	
		printf("\n\t Using device:\t\t\t%s\n", devProp.name);
		printf("\n\t Max grid size:\t\t\t%d\n", devProp.maxGridSize[0]);
		maxgrid_x = devProp.maxGridSize[0];
	}
		checkCudaErrors(cudaSetDevice(device));
	
//------------ memory setup -------------------------------------
	
	// if set before set device getting errors - invalid handle
	GpuTimer timer;  

	float2 *d_spectra;
	float  *d_real, *d_img, *d_coeff;

	float fir_time = 0.0f;
	float fft_time = 0.0f;
	float mem_time_in = 0.0f;
	float mem_time_out = 0.0f;

	//malloc
	printf("\nDevice memory allocation...\t\t");
	timer.Start();
	checkCudaErrors(cudaMalloc((void **) &d_spectra, sizeof(float2)*(filesize)));
	checkCudaErrors(cudaMalloc((void **) &d_coeff,   sizeof(float)*nChannels*nTaps));
	checkCudaErrors(cudaMalloc((void **) &d_real,    sizeof(float)*filesize));
	checkCudaErrors(cudaMalloc((void **) &d_img,     sizeof(float)*filesize));
	timer.Stop();
	printf("done in %g ms.", timer.Elapsed());

	printf("\n\t\td_spectra using filesize: \t%g MB.", sizeof(float2)*(filesize)/1024.0/1024);
	printf("\n\t\td_coeff using filesize: \t%g MB.", sizeof(float)*filesize/1024.0/1024);
	printf("\n\t\td_real using filesize: \t\t%g MB.", sizeof(float)*filesize/1024.0/1024);
	printf("\n\t\td_img using filesize: \t\t%g MB.", sizeof(float)*filesize/1024.0/1024);
	printf("\n\t\t----------------------\t\t-----------");
	printf("\n\t\tTotal: \t\t\t\t%g MB.\n\n",sizeof(float)*(5.0*filesize)/1024/1024);
	
	// set to 0.0
	printf("\nDevice memset...\t\t\t");
		timer.Start();
			checkCudaErrors(cudaMemset(d_spectra, 0.0, sizeof(float2)*(filesize)));
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
	mem_time_in=timer.Elapsed();
	printf("done in %g ms.\n", mem_time_in);

	//--------------- Fir ----------------------------
		
	int run_blocks = nBlocks - nTaps + 1; 
	int grid_x;
	int n_cycle = run_blocks/maxgrid_x + 1;
	printf("n_cycle : %d \n", n_cycle);

	printf("\n\t\t------------ Kernel run-----------------");
	for (int i = 0; i < n_cycle; i++){
		
		if (maxgrid_x < run_blocks){
									grid_x = maxgrid_x;
									run_blocks = run_blocks - maxgrid_x;
															
		} else grid_x = run_blocks;
		
		dim3 gridSize(grid_x, blocks_y, 1);
		dim3 blockSize(nChannels/gridSize.y, 1, 1); 
	
		timer.Start();
			Fir<<<gridSize, blockSize>>>(d_real + i*(maxgrid_x)*nChannels, d_img + i*(maxgrid_x)*nChannels, d_coeff, nTaps, nChannels, d_spectra + i*maxgrid_x*nChannels);
		timer.Stop();
		fir_time+=timer.Elapsed();
	
		//----- error check -----
			checkCudaErrors(cudaGetLastError());
			checkCudaErrors(cudaDeviceSynchronize());
		//-----------------------
	
		printf("\nFir kernel %d\n", i);
		printf("\n\n blocks \t time \t\t threads \t bandwidth \t flops");
		printf("\n%d \t\t %lf \t %i \t\t %g \t %g\n",grid_x, fir_time, nChannels/blocks_y, 
													(3*nChannels*nTaps*sizeof(float)+nChannels*2*sizeof(float))*(nBlocks-nTaps+1)*1000/fir_time, 
													(4*nTaps)*nChannels*(nBlocks-nTaps+1)*1000.0/fir_time);
}
//--------------- cuFFT ----------------------------
/*
	//Create fft Plan
	cufftHandle plan;
	cufftPlan1d(&plan, nChannels, CUFFT_C2C, nBlocks);
	
	//execute plan and copy back to host
	printf("\n\ncuFFT run..\t\t");
		timer.Start();
			cufftExecC2C(plan, (cufftComplex *)d_spectra, (cufftComplex *)d_spectra, CUFFT_FORWARD);
		timer.Stop();
		fft_time = timer.Elapsed();
	printf("done in %g ms.\n\n", timer.Elapsed());
	
	//Destroy the cuFFT plan
	cufftDestroy(plan);
*/


//--------------- copy data back ----------------------------
	printf("Copy data from device to host \t");
	timer.Start();
		checkCudaErrors(cudaMemcpy(spectra,d_spectra,(filesize)*sizeof(float2), cudaMemcpyDeviceToHost));	
	timer.Stop();
	printf("done in %g ms.\n", timer.Elapsed());
	mem_time_out=timer.Elapsed();

printf("\nTotal execution time %g ms.\n", mem_time_in + fir_time + mem_time_out);

//---------------- write to file process ----------------------

	char str[200];
	sprintf(str,"GPU-polyphase.dat");
	
		printf("\n Write results into file...\t");
		if (WRITE) save_time(str, nBlocks-nTaps+1, fir_time, fft_time, mem_time_in, mem_time_out, nChannels, nTaps);
		printf("\t done.\n-------------------------------------\n");

//--------------- clean-up process ----------------------------
	
	checkCudaErrors(cudaFree(d_spectra));
	checkCudaErrors(cudaFree(d_real));
	checkCudaErrors(cudaFree(d_img));
	checkCudaErrors(cudaFree(d_coeff));
	// device reset is in the polyphase-analyze.c
	
}
