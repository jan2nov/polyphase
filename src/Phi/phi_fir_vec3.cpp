#include <mkl.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>


double FIR_check_uni(float *input_data, float *spectra_GPU, float *coeff, unsigned int nTaps, unsigned int nChannels, unsigned int nBlocks){
	nChannels=nChannels*2;
	float ftemp=0;
	float etemp=0;
	double error=0;
	int count=0;
	//cout << "nBlocks: " << nBlocks-nTaps+1 << endl;
	for(int bl=0;bl<nBlocks-nTaps+1;bl++){
		for(int c=0;c<nChannels;c++){
			ftemp=0;
			for(int t=0;t<nTaps;t++){
				ftemp+=coeff[t*nChannels + c]*input_data[bl*nChannels + t*nChannels + c];
			}//nTaps
			count++;
			etemp=(ftemp-spectra_GPU[bl*nChannels + c])*(ftemp-spectra_GPU[bl*nChannels + c]);
			error+=etemp;
		}//nChannels
	} //nBlocks

	return(error);
}


//lapack_complex_float==MKL_Complex8 = (real,float)
//http://software.intel.com/sites/products/documentation/hpc/mkl/lin/MKL_UG_coding_calls/Using_Complex_Types_in_C_C.htm

void Phi_FFT_separated_single(MKL_Complex8 *input, int nSpectra, int nChannels){
	int nThreads = omp_get_num_threads();
	int nTransformations=nSpectra/nThreads;
	int th=0;
	
	DFTI_DESCRIPTOR_HANDLE handle = 0;
	DftiCreateDescriptor(&handle, DFTI_SINGLE, DFTI_COMPLEX,1,(MKL_LONG) nChannels);
	DftiSetValue(handle, DFTI_PLACEMENT, DFTI_INPLACE); //I will over write input array
	DftiCommitDescriptor(handle);	
	//schedule(static,block_step)
	#pragma omp parallel for shared(nTransformations,nChannels,nThreads,input) private(temp)
	for (th = 0; th < nSpectra; th++) {
		DftiComputeForward(handle, &input[th*nChannels]);
	}
	DftiFreeDescriptor(&handle);
}


void FIR_XEON_PHI( int nTaps, int nChannels, int nSpectra){
	//********************************* Allocations **********************************
	double FIR_time,FIR_start,FFT_time,FFT_start;
	double Total_time=0;
	long int bl,t;

	int nBlocks=nSpectra+nTaps-1;

	unsigned int data_size=nBlocks*nChannels*2;
	unsigned int spectra_size=nSpectra*nChannels*2;
	unsigned int coeff_size=nTaps*nChannels*2;	
	
	static float  *input_data;
	static float  *spectra;
	static float  *coeff;
	
	input_data = (float*)_mm_malloc( data_size*sizeof(float) ,64);
	spectra = (float*)_mm_malloc( spectra_size*sizeof(float) ,64);
	coeff = (float*)_mm_malloc( coeff_size * sizeof(float) ,64);
	//********************************* Allocations **********************************
	
	
	// initialize
	#pragma omp parallel for
	for(t=0;t<data_size;t++) {
		input_data[t]=(int) t/(nChannels*2);
	}
	#pragma omp parallel for
	for(t=0;t<spectra_size;t++) {
		spectra[t]=0;
	}
	#pragma omp parallel for
	for(t=0;t<coeff_size;t++) {
		coeff[t]=1.0;
	}
		
	FIR_start=omp_get_wtime();
	
	
	#pragma omp parallel for private(t,bl) shared(input_data,spectra,coeff)
	for(bl=0; bl<(nBlocks-nTaps+1); bl++) {
		//block=bl*nChannels*2;
		for(t=0;t<nTaps;t++){
			// using Cilk vector notation
			#pragma vector aligned (input_data,spectra,coeff)
			spectra[bl*nChannels*2:nChannels*2]=coeff[t*nChannels*2:nChannels*2]*input_data[(t+bl)*nChannels*2:nChannels*2] + spectra[bl*nChannels*2:nChannels*2];
		}
	}
	FIR_time = omp_get_wtime() - FIR_start;
	
	//printf("Error: %f \n", FIR_check_uni(input_data,spectra,coeff,nTaps,nChannels,nBlocks));
	
	FFT_start=omp_get_wtime();
	Phi_FFT_separated_single((MKL_Complex8 *) spectra,nSpectra,nChannels);
	FFT_time = omp_get_wtime() - FFT_start;
	

	double Gflop=2.0*nTaps*data_size;
	double GFlops=Gflop/FIR_time;
	double Bandwidth=(8.0*data_size*nTaps+4.0*data_size)/FIR_time;
	
	double FFTGFlops=nSpectra*(2.5*nChannels*log2(1.0*nChannels))/FFT_time;
	
	Total_time=FIR_time+FFT_time;
	printf("%d %f %f %f %f %f %f\n",nSpectra,FIR_time,Total_time,Bandwidth,GFlops,FFTGFlops,FFT_time);
	
	
	
	_mm_free(input_data);
	_mm_free(spectra);
	_mm_free(coeff);
}



int main(){
	int nTaps=8;
	int nChannels=512;
	int nSpectra=120000;

	double FIR_time,Total_time;
	long int bl,t,i;
	int data_start;
	
	//omp_set_num_threads(120);
	//kmp_set_defaults("KMP_AFFINITY=compact");
	
	printf("Starting Calculation\n");
	//printf("nSpectra=%d ;nBlocks=%d\n",nSpectra,nBlocks);
	//printf("Data size: %d ;Spectra size: %d \n",data_size,spectra_size);
	
	nSpectra=15000;
	FIR_XEON_PHI(nTaps,nChannels,nSpectra);
	
	nSpectra=15000;
	FIR_XEON_PHI(nTaps,nChannels,nSpectra);
	
	nSpectra=30000;
	FIR_XEON_PHI(nTaps,nChannels,nSpectra);
	
	nSpectra=60000;
	FIR_XEON_PHI(nTaps,nChannels,nSpectra);
	
	nSpectra=120000;
	FIR_XEON_PHI(nTaps,nChannels,nSpectra);
	
	nSpectra=240000;
	FIR_XEON_PHI(nTaps,nChannels,nSpectra);	
	
	nSpectra=480000;
	FIR_XEON_PHI(nTaps,nChannels,nSpectra);	
	
	//printf("duration=%f result=%f\n",duration,spectra[0]);
	//printf("Running %d openmp threads \n", omp_get_max_threads());
	//printf("SP Gflops = %f\n",Gflops);
		
	return (0);
}