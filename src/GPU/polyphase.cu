#include <cuda.h>
#include <stdio.h>
#include <cufft.h>
#include <string.h>
#include "utils_cuda.h"
#include "utils_file.h"
//#include <cutil_inline.h>
#include "data.h"
#include "timer.h"

#define NUMBER 3
#define NTIMES 1
y
__constant__ float C_coeff[4096];

__global__ void Fir( const float *d_signal_real, const float *d_signal_img, const float *coeff, const int nTaps, const int nChannels, float2 *spectra)
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

__global__ void Fir_SpB(float2* __restrict__  d_data, float* __restrict__ d_coeff, int nTaps, int nChannels, int yshift, float2* __restrict__ d_spectra) {
	int t = 0;
	int bl= blockIdx.x*nChannels;
	int ypos = blockDim.x*blockIdx.y + yshift;
	float2 ftemp1 = make_float2(0.0,0.0);

	for(t=ypos + threadIdx.x;t<(nTaps)*nChannels;t+=nChannels){
		ftemp1.x  += d_coeff[t]*d_data[bl+t].x;
		ftemp1.y  += d_coeff[t]*d_data[bl+t].y;
	}

	t=bl + ypos + threadIdx.x;
	d_spectra[t]=ftemp1;
	return;
}
__global__ void Fir_prefetch(float2* __restrict__  d_data, float* __restrict__ d_coeff, int nTaps, int nChannels, int yshift, float2* __restrict__ d_spectra) {
	int t = 0;
	int bl= blockIdx.x*nChannels;
	int ypos = blockDim.x*blockIdx.y + yshift;
	float2 ftemp1 = make_float2(0.0,0.0);
	float2 ftemp2 = make_float2(0.0,0.0);
	float2 sum = make_float2(0.0,0.0);
	//ftemp1.x=0.0f;ftemp1.y=0.0f;

	ftemp1.x  = d_coeff[ypos+threadIdx.x]*d_data[bl+ypos+threadIdx.x].x;
	ftemp1.y  = d_coeff[ypos+threadIdx.x]*d_data[bl+ypos+threadIdx.x].y;
	for(t=ypos + threadIdx.x;t<(nTaps-1)*nChannels;t+=nChannels){
		ftemp2.x  = d_coeff[t+nChannels]*d_data[bl+t+nChannels].x;
		ftemp2.y  = d_coeff[t+nChannels]*d_data[bl+t+nChannels].y;
		sum.x += ftemp1.x;
		sum.y += ftemp1.y;
		ftemp1 = ftemp2;
	}
	sum.x += ftemp1.x;
	sum.y += ftemp1.y;

	t=bl + ypos + threadIdx.x;
	d_spectra[t]=sum;
	return;
}

__global__ void Fir_simple(float2* d_data, float* d_coeff, int nTaps, int nChannels, int yshift, float2* d_spectra) {
	int t = 0;
	int bl= 5*blockIdx.x*nChannels;
	int ypos = blockDim.x*blockIdx.y + yshift;
	float2 ftemp1 = make_float2(0.0,0.0);
	float2 ftemp2 = make_float2(0.0,0.0);
	float2 ftemp3 = make_float2(0.0,0.0);
	float2 ftemp4 = make_float2(0.0,0.0);
	float2 ftemp5 = make_float2(0.0,0.0);
	float temp;

	for(t=ypos + threadIdx.x;t<(nTaps)*nChannels;t+=nChannels){
	  temp = d_coeff[t]; 
		ftemp1.x += temp*d_data[bl+t].x;
		ftemp1.y += temp*d_data[bl+t].y;
		ftemp2.x += temp*d_data[bl+nChannels+t].x;
		ftemp2.y += temp*d_data[bl+nChannels+t].y;
		ftemp3.x += temp*d_data[bl+2*nChannels+t].x;
		ftemp3.y += temp*d_data[bl+2*nChannels+t].y;
		ftemp4.x += temp*d_data[bl+3*nChannels+t].x;
		ftemp4.y += temp*d_data[bl+3*nChannels+t].y;
		ftemp5.x += temp*d_data[bl+4*nChannels+t].x;
		ftemp5.y += temp*d_data[bl+4*nChannels+t].y;
	}

	t=bl + ypos + threadIdx.x;
	d_spectra[t]=ftemp1;
	d_spectra[t+nChannels]=ftemp2;
	d_spectra[t+2*nChannels]=ftemp3;
	d_spectra[t+3*nChannels]=ftemp4;
	d_spectra[t+4*nChannels]=ftemp5;
	return;
}

__global__ void Fir_restrict(float2* __restrict__ const d_data, float* __restrict__ const d_coeff, int nTaps, int nChannels, int yshift, float2* __restrict__ d_spectra) {
	int t = 0;
	int bl= 4*blockIdx.x*nChannels;
	int ypos = blockDim.x*blockIdx.y + yshift;
	float2 ftemp1 = make_float2(0.0,0.0);
	float2 ftemp2 = make_float2(0.0,0.0);
	float2 ftemp3 = make_float2(0.0,0.0);
	float2 ftemp4 = make_float2(0.0,0.0);
	float temp;

	for(t=ypos + threadIdx.x;t<(nTaps)*nChannels;t+=nChannels){
	  temp = d_coeff[t]; 
		ftemp1.x += temp*d_data[bl+t].x;
		ftemp1.y += temp*d_data[bl+t].y;
		ftemp2.x += temp*d_data[bl+nChannels+t].x;
		ftemp2.y += temp*d_data[bl+nChannels+t].y;
		ftemp3.x += temp*d_data[bl+2*nChannels+t].x;
		ftemp3.y += temp*d_data[bl+2*nChannels+t].y;
		ftemp4.x += temp*d_data[bl+3*nChannels+t].x;
		ftemp4.y += temp*d_data[bl+3*nChannels+t].y;
	}

	t=bl + ypos + threadIdx.x;
	d_spectra[t]=ftemp1;
	d_spectra[t+nChannels]=ftemp2;
	d_spectra[t+2*nChannels]=ftemp3;
	d_spectra[t+3*nChannels]=ftemp4;
	return;
}

__global__ void Fir_ldg(float2 const* __restrict__ d_data, float const* __restrict__ d_coeff, int nTaps, int nChannels, int yshift, float2* d_spectra) {
	int t = 0;
	int bl= NUMBER*blockIdx.x*nChannels;
	int ypos = blockDim.x*blockIdx.y + yshift;
	float temp;
	float2 ftemp[NUMBER];
	
	for (int i = 0; i<NUMBER;i++){
	  ftemp[i].x = 0.0f;
	  ftemp[i].y = 0.0f;
	}

	for(t=yshift + threadIdx.x+ypos; t<(nTaps)*nChannels; t+=nChannels){
	  //temp = __ldg(&d_coeff[t]); 
  	  temp = d_coeff[t]; 
	  for(int i=0;i<NUMBER;i++){
	    //ftemp[i].x += temp*__ldg(&d_data[bl+i*nChannels + t].x);
	    ftemp[i].x += temp*(d_data[bl+i*nChannels + t].x);
	    //ftemp[i].y += temp*__ldg(&d_data[bl+i*nChannels + t].y);
	    ftemp[i].y += temp*(d_data[bl+i*nChannels + t].y);
	  }
	}

	t=bl + yshift + threadIdx.x+ypos;
	for (int i=0;i<NUMBER;i++){
	  d_spectra[t + i*nChannels]=ftemp[i];
	}

	return;
}
/*
__global__ void Fir_ldg2(float2* d_data, float* d_coeff, int nTaps, int nChannels, int yshift, float2* d_spectra) {
	int t = 0;
	int bl= 5*blockIdx.x*nChannels;
	//int pos = 2*gridDim.x;
	float2 ftemp1 = make_float2(0.0,0.0);
	float2 ftemp2 = make_float2(0.0,0.0);
	float2 ftemp3 = make_float2(0.0,0.0);
	float2 ftemp4 = make_float2(0.0,0.0);
	float2 ftemp5 = make_float2(0.0,0.0);
	float temp;

	for(t=yshift + threadIdx.x;t<(nTaps)*nChannels;t+=nChannels){
	  temp = __ldg(&d_coeff[t]); 
	  ftemp1.x += temp*__ldg(&d_data[bl+t].x);
	  ftemp1.y += temp*__ldg(&d_data[bl+t].y);
	  ftemp2.x += temp*__ldg(&d_data[bl+nChannels+t].x);
	  ftemp2.y += temp*__ldg(&d_data[bl+nChannels+t].y);
	  ftemp3.x += temp*__ldg(&d_data[bl+2*nChannels+t].x);
	  ftemp3.y += temp*__ldg(&d_data[bl+2*nChannels+t].y);
	  ftemp4.x += temp*__ldg(&d_data[bl+3*nChannels+t].x);
	  ftemp4.y += temp*__ldg(&d_data[bl+3*nChannels+t].y);
	  ftemp5.x += temp*__ldg(&d_data[bl+4*nChannels+t].x);
	  ftemp5.y += temp*__ldg(&d_data[bl+4*nChannels+t].y);
	}

	t=bl + yshift + threadIdx.x;
	d_spectra[t]=ftemp1;
	d_spectra[t+nChannels]=ftemp2;
	d_spectra[t+2*nChannels]=ftemp3;
	d_spectra[t+3*nChannels]=ftemp4;
	d_spectra[t+4*nChannels]=ftemp5;
	return;
}
*/
__global__ void Fir_fmaf(float2* __restrict__ const d_data, float* __restrict__ const d_coeff, int nTaps, int nChannels, int yshift, float2* __restrict__ d_spectra) {
	int t = 0;
	int bl= NUMBER*blockIdx.x*nChannels;
	int ypos = blockDim.x*blockIdx.y + yshift;
	float2 ftemp[NUMBER];
	float temp;
	
	for (int i = 0; i<NUMBER;i++){
	  ftemp[i].x = 0.0f;
	  ftemp[i].y = 0.0f;
	}
		
	for(t=ypos + threadIdx.x;t<(nTaps)*nChannels;t+=nChannels){
	  temp = d_coeff[t]; 
	    for(int i=0;i<NUMBER;i++){
		ftemp[i].x = __fmaf_rn(temp,d_data[bl+i*nChannels + t].x, ftemp[i].x);
		ftemp[i].y = __fmaf_rn(temp,d_data[bl+i*nChannels + t].y, ftemp[i].y);
		}
	}

	t=bl + yshift + threadIdx.x;
	for (int i=0;i<NUMBER;i++){
	  d_spectra[t + i*nChannels]=ftemp[i];
	}
	
	return;
}
/*
__global__ void Fir_fmaf_ldg(float2* d_data, float* d_coeff, int nTaps, int nChannels, int yshift, float2* d_spectra) {
	int t = 0;
	int bl= 4*blockIdx.x*nChannels;
	int ypos = blockDim.x*blockIdx.y + yshift;
	float2 ftemp1 = make_float2(0.0,0.0);
	float2 ftemp2 = make_float2(0.0,0.0);
	float2 ftemp3 = make_float2(0.0,0.0);
	float2 ftemp4 = make_float2(0.0,0.0);
	float temp;

	for(t=ypos + threadIdx.x;t<(nTaps)*nChannels;t+=nChannels){
	  temp = __ldg(&d_coeff[t]); 
		ftemp1.x = __fmaf_rn(temp,__ldg(&d_data[bl+t].x),ftemp1.x);
		ftemp1.y = __fmaf_rn(temp,__ldg(&d_data[bl+t].y),ftemp1.y);
		ftemp2.x = __fmaf_rn(temp,__ldg(&d_data[bl+nChannels+t].x),ftemp2.x);
		ftemp2.y = __fmaf_rn(temp,__ldg(&d_data[bl+nChannels+t].y),ftemp2.y);
		ftemp3.x = __fmaf_rn(temp,__ldg(&d_data[bl+2*nChannels+t].x),ftemp3.x);
		ftemp3.y = __fmaf_rn(temp,__ldg(&d_data[bl+2*nChannels+t].y),ftemp3.y);
		ftemp4.x = __fmaf_rn(temp,__ldg(&d_data[bl+3*nChannels+t].x),ftemp4.x);
		ftemp4.y = __fmaf_rn(temp,__ldg(&d_data[bl+3*nChannels+t].y),ftemp4.y);
	}

	t=bl + ypos + threadIdx.x;
	d_spectra[t]=ftemp1;
	d_spectra[t+nChannels]=ftemp2;
	d_spectra[t+2*nChannels]=ftemp3;
	d_spectra[t+3*nChannels]=ftemp4;
	return;
}
*/
void gpu_code(  float2 *data_in, 
		float2 *spectra, 
		float *coeff,
		const int nChannels,
		unsigned int nBlocks, 
		unsigned int filesize,
		const int nThreads,
		int blocks_y,
		const int nTaps){

 bool WRITE=true;					
//------------ initialize card -----------

  int devCount, device;
  int maxgrid_x = 0;
  
  checkCudaErrors(cudaGetDeviceCount(&devCount));
  printf("\n\t\t-------------- GPU part -----------------");
  printf("\nThere are %d devices.", devCount);

	for (int i = 0; i < devCount; i++){
		cudaDeviceProp devProp;
		device = 0;
		checkCudaErrors(cudaGetDeviceProperties(&devProp,device));	
		printf("\n\t Using device:\t\t\t%s\n", devProp.name);
		printf("\n\t Max grid size:\t\t\t%d\n", devProp.maxGridSize[0]);
		maxgrid_x = devProp.maxGridSize[0]-3;
	}
		checkCudaErrors(cudaSetDevice(device));
	
//------------ memory setup -------------------------------------
	
	// if set before set device getting errors - invalid handle
	GpuTimer timer;  

	float2 *d_spectra, *d_data_in;
	float  *d_coeff;

	float fir_time = 0.0f;
	float fft_time = 0.0f;
	float mem_time_in = 0.0f;
	float mem_time_out = 0.0f;
	unsigned int orig_filesize = filesize;
	unsigned int orig_nBlocks = nBlocks;
	int nSteps = 1;
	size_t mem_free_0=0, mem_tot_0=0, mem_free_1=0;
	
	checkCudaErrors(cudaMemGetInfo(&mem_free_0, & mem_tot_0));
	printf("\nFree memory: %zu \t Total Memory: %zu\n", mem_free_0/1024/1024, mem_tot_0/1024/1024);

	if (filesize > (60000+nTaps-1)*nChannels) {
	   filesize = (60000+nTaps-1)*nChannels;
	   nBlocks = filesize/nChannels;
	   nSteps = (orig_nBlocks - nTaps + 1)/(nBlocks - nTaps + 1);
	}

	//malloc
	printf("\nDevice memory allocation...: \t\t");
	timer.Start();
	checkCudaErrors(cudaMalloc((void **) &d_spectra, sizeof(float2)*(orig_filesize)));
	checkCudaErrors(cudaMalloc((void **) &d_coeff,   sizeof(float)*nChannels*nTaps));
	checkCudaErrors(cudaMalloc((void **) &d_data_in,    sizeof(float2)*filesize));
	timer.Stop();
	printf("done in %g ms.", timer.Elapsed());
	checkCudaErrors(cudaMemGetInfo(&mem_free_1, &mem_tot_0));

	printf("\n\t\tFilesize orig: \t%u spektra.", orig_filesize/nChannels);
	printf("\n\t\td_spectra using filesize: \t%g MB.", sizeof(float2)*(orig_filesize)/1024.0/1024);
	printf("\n\t\td_coeff using filesize: \t%g MB.", sizeof(float)*nTaps*nChannels/1024.0/1024);
	printf("\n\t\td_real using filesize: \t\t%g MB.", sizeof(float)*filesize/1024.0/1024);
	printf("\n\t\td_img using filesize: \t\t%g MB.", sizeof(float)*filesize/1024.0/1024);
	printf("\n\t\t----------------------\t\t-----------");
	printf("\n\t\tTotal: \t\t\t\t%g MB.",sizeof(float)*(4.0*filesize)/1024/1024+sizeof(float)*nTaps*nChannels/1024.0/1024);
	printf("\n\t\tAllocated: \t\t\t%zu MB.\n\n", (mem_free_0-mem_free_1)/1024/1024);
	
	// set to 0.0
	printf("\nDevice memset...\t\t\t");
		timer.Start();
			checkCudaErrors(cudaMemset(d_spectra, 0.0, sizeof(float2)*(orig_filesize)));
		timer.Stop();
	printf("done in %g ms.", timer.Elapsed());

	// copy data to device
	for (int j=0; j<nSteps*(nBlocks-nTaps+1); j+=nBlocks-nTaps+1){
	printf("\n-------------- run %i------------\n", j);
	printf("\nCopy data from host to device...\t");
		timer.Start();
			checkCudaErrors(cudaMemcpy(d_coeff, coeff, nChannels*nTaps*sizeof(float), cudaMemcpyHostToDevice));
			//checkCudaErrors(cudaMemcpyToSymbol(C_coeff,  coeff,   sizeof(float)*nChannels*nTaps));
			checkCudaErrors(cudaMemcpy(d_data_in, data_in+j*nChannels, filesize*sizeof(float2), cudaMemcpyHostToDevice));
		timer.Stop();
	mem_time_in+=timer.Elapsed();
	printf("done in %g ms.\n", mem_time_in);

	printf("\n\t\t------------ Kernel run-----------------\n");		
	int run_blocks = (int)(nBlocks - nTaps + 1); 
	int grid_x;
	int n_cycle = run_blocks/maxgrid_x + 1;
	printf("n_cycle : %d \t nTaps: %d\n", n_cycle, nTaps);
	
	cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
 for (int k = 0; k<NTIMES ; k++){
	for (int i = 0; i < n_cycle; i++){
		
		if (maxgrid_x < run_blocks){
									grid_x = maxgrid_x;
									run_blocks = run_blocks - maxgrid_x;
		} else grid_x = run_blocks;
		
		dim3 gridSize(grid_x/NUMBER, blocks_y, 1);
		dim3 blockSize(nThreads/gridSize.y, 1, 1); 
	
		timer.Start();
		int nKernels=1;//(int) nChannels/blockSize.x;
		int remainder=0;//nChannels-nKernels*blockSize.x;
		for (int nutak=0;nutak<nKernels;nutak++){	
		  Fir_ldg<<<gridSize, blockSize>>>((float2*) d_data_in + i*(maxgrid_x)*nChannels, d_coeff, nTaps, nChannels, nutak*blockSize.x, (float2*) d_spectra + i*maxgrid_x*nChannels+j*nChannels);
		}
		if (remainder>0){
			int old_blockSize = blockSize.x;
			blockSize.x=remainder;
			Fir_ldg<<<gridSize, blockSize>>>((float2*) d_data_in + i*(maxgrid_x)*nChannels, d_coeff, nTaps, nChannels, nKernels*old_blockSize, (float2*) d_spectra + i*maxgrid_x*nChannels+j*nChannels);
		}
		timer.Stop();
	 	fir_time+=timer.Elapsed();
		/*timer.Start();
			Fir<<<gridSize, blockSize>>>(d_real + i*(maxgrid_x)*nChannels, d_img + i*(maxgrid_x)*nChannels, d_coeff, nTaps, nChannels, d_spectra + i*maxgrid_x*nChannels);
		timer.Stop();
		fir_time+=timer.Elapsed();
		*/
		//----- error check -----
			checkCudaErrors(cudaGetLastError());
			//checkCudaErrors(cudaDeviceSynchronize());
		//-----------------------

		printf("\nFir kernel %d\n", i);
		printf("\n\n blocks \t time \t\t threads \t bandwidth \t flops");
		printf("\n%d \t\t %lf \t %i \t\t %g \t %g\n",grid_x, fir_time, nThreads/blocks_y, 
													(3.0*nChannels*nTaps*sizeof(float)+nChannels*2.0*sizeof(float))*(orig_nBlocks-nTaps+1)*1000.0*NTIMES/fir_time, 
													(4*nTaps)*nChannels*(orig_nBlocks-nTaps+1)*1000.0*NTIMES/fir_time);
}
}
}
//--------------- cuFFT ----------------------------

checkCudaErrors(cudaFree(d_data_in));
checkCudaErrors(cudaFree(d_coeff));

	//Create fft Plan
	cufftHandle plan;
	cufftResult error;
	error = cufftPlan1d(&plan, nChannels, CUFFT_C2C, orig_nBlocks);
	if (CUFFT_SUCCESS != error){
		printf("CUFFT error: %d", error);
	}
	
	//execute plan and copy back to host
	printf("\n\ncuFFT run..\t\t");
		timer.Start();
			cufftExecC2C(plan, (cufftComplex *)d_spectra, (cufftComplex *)d_spectra, CUFFT_FORWARD);
		timer.Stop();
		fft_time = timer.Elapsed();
	printf("done in %g ms.\n\n", timer.Elapsed());
	
	//Destroy the cuFFT plan
	cufftDestroy(plan);



//--------------- copy data back ----------------------------
	printf("Copy data from device to host \t");
	timer.Start();
		checkCudaErrors(cudaMemcpy(spectra,d_spectra,(orig_filesize)*sizeof(float2), cudaMemcpyDeviceToHost));	
	timer.Stop();
	printf("done in %g ms.\n", timer.Elapsed());
	mem_time_out=timer.Elapsed();

printf("\nTotal execution time %g ms.\n", mem_time_in + fir_time + mem_time_out);

//---------------- write to file process ----------------------

	char str[200];
	sprintf(str,"GPU-polyphase.dat");
	
		printf("\n Write results into file...\t");
		if (WRITE) save_time(str, (orig_nBlocks-nTaps+1)*NTIMES, fir_time, fft_time, mem_time_in, mem_time_out, nChannels, nTaps, 0, nThreads);
		printf("\t done.\n-------------------------------------\n");

//--------------- clean-up process ----------------------------
	
	checkCudaErrors(cudaFree(d_spectra));
	// device reset is in the polyphase-analyze.c
	
}
