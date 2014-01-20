#include <omp.h>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <time.h>
#include <fftw3.h>
#include <string.h>
#include <xmmintrin.h>

#include <sys/time.h>
#include <unistd.h>

#include "Filter_window.h"
#include "Filters.h"

using namespace std;

bool check_errors=0;

void FFT(float *spectra_real, float *spectra_img, unsigned int nTaps, unsigned int nChannels, unsigned int nBlocks){
	unsigned int f,bl,th_id,nthreads,block_step;
	fftwf_complex *in, *out;
	fftwf_plan p;
	
	#pragma omp parallel private(in,out,p)
	{
		th_id = omp_get_thread_num();
		nthreads = omp_get_num_threads();
		block_step=nBlocks/nthreads;
		in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChannels);
		out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChannels);
		#pragma omp critical (p)
		p = fftwf_plan_dft_1d(nChannels, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		if (!p){
			exit(101);
		}
		
		#pragma omp for schedule(static,block_step)
		for(bl=0;bl<nBlocks;bl++){
			// I need better memory arrangement
			for(f=0;f<nChannels;f++){
				in[f][0]=spectra_real[bl*nChannels+f];
				in[f][1]=spectra_img[bl*nChannels+f];
			}
			fftwf_execute(p);
			for(f=0;f<nChannels;f++){
				spectra_real[bl*nChannels+f]=out[f][0];
				spectra_img[bl*nChannels+f]=out[f][1];
			}
		}
		
		fftwf_destroy_plan(p);
		fftwf_free(in);
		fftwf_free(out);
	}
}


void FFT_mod(float *spectra, unsigned int nTaps, unsigned int nChannels, unsigned int nBlocks){
	unsigned int f,bl,th_id,nthreads,block_step;
	fftwf_complex *in, *out;
	fftwf_plan p;
	
	#pragma omp parallel private(in,out,p)
	{
		th_id = omp_get_thread_num();
		nthreads = omp_get_num_threads();
		block_step=nBlocks/nthreads;
		in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChannels);
		out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChannels);
		#pragma omp critical (p)
		p = fftwf_plan_dft_1d(nChannels, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		if (!p){
			exit(101);
		}
		
		#pragma omp for schedule(static,block_step)
		for(bl=0;bl<nBlocks;bl++){
			// I need better memory arrangement
			for(f=0;f<nChannels;f++){
				in[f][0]=spectra[bl*nChannels*2 + 2*f];
				in[f][1]=spectra[bl*nChannels*2 + 2*f + 1];
			}
			fftwf_execute(p);
			for(f=0;f<nChannels;f++){
				spectra[bl*nChannels*2 + 2*f]=out[f][0];
				spectra[bl*nChannels*2 + 2*f + 1]=out[f][1];
			}
		}
		
		fftwf_destroy_plan(p);
		fftwf_free(in);
		fftwf_free(out);
	}
}

//************************************************************************************************************
//*                 DATA CHECKING
//************************************************************************************************************
double FIR_check_sep(float *real, float *img, float *spectra_real_GPU, float *spectra_img_GPU, float *coeff, unsigned int nTaps, unsigned int nChannels, unsigned int nBlocks){
	float ftemp_real=0,ftemp_img=0;
	double error_real=0,error_img=0;
	for(int bl=0;bl<nBlocks-nTaps+1;bl++){
		for(int c=0;c<nChannels;c++){
			ftemp_real=0;ftemp_img=0;
			for(int t=0;t<nTaps;t++){
				ftemp_real+=coeff[t*nChannels + c]*real[bl*nChannels + t*nChannels + c];
				ftemp_img+=coeff[t*nChannels + c]*img[bl*nChannels + t*nChannels + c];
			}//nTaps
			error_real+=(ftemp_real-spectra_real_GPU[bl*nChannels + c])*(ftemp_real-spectra_real_GPU[bl*nChannels + c]);
			error_img+=(ftemp_img-spectra_img_GPU[bl*nChannels + c])*(ftemp_img-spectra_img_GPU[bl*nChannels + c]);
		}//nChannels
	} //nBlocks

	return(error_real+error_img);
}

double FIR_check_uni(float *input_data, float *spectra_GPU, float *coeff, unsigned int nTaps, unsigned int nChannels, unsigned int nBlocks){
	nChannels=nChannels*2;
	float ftemp=0;
	float etemp=0;
	double error=0;
	unsigned int count=0;
	for(int bl=0;bl<nBlocks-nTaps+1;bl++){
		for(int c=0;c<nChannels;c++){
			ftemp=0;
			for(int t=0;t<nTaps;t++){
				ftemp+=coeff[t*nChannels + c]*input_data[bl*nChannels + t*nChannels + c];
			}//nTaps
			etemp=(ftemp-spectra_GPU[bl*nChannels + c])*(ftemp-spectra_GPU[bl*nChannels + c]);
			error+=etemp;
		}//nChannels
	} //nBlocks

	return(error);
}

//************************************************************************************************************
//*                 DATA CHECKING
//************************************************************************************************************


void Test_intrinsic_Mk3(long int DATA_SIZE, int nTaps, int nChannels, double *mtime_Decoupling, double *mtime_FIR, double *mtime_FFT, double *mtime_Coupling, double *mtime){
	long int data_size=0,separated_data_size=0,f=0;
	double start, end, Decoupling_start, Decoupling_end, FIR_start, FIR_end, FFT_start, FFT_end, Coupling_start, Coupling_end;
	float *input_data;
	float *coeff_mod;
	float *spectra;
	int nBlocks;
	srand (time(NULL));
	
	data_size=DATA_SIZE*2;
	input_data=new float[data_size];
	spectra=new float[data_size];
	
	for(f=0;f<data_size;f++){
		input_data[f]=(float) (rand() % 10000 + 1)/1000.0;
	}
	
	separated_data_size=data_size/2;
	nBlocks=separated_data_size/nChannels; // determining number of blocks
	
	//***************************** PERFORMANCE MEASUREMENTS *******************************
	start=omp_get_wtime();
	
	//********************************* Window *********************************************
	coeff_mod = new float [nChannels * nTaps * 2];
	Window_coefficients_mod(nTaps,nChannels,1,coeff_mod);
	//********************************* Window *********************************************
	
	//********************************* FIR *********************************************
	FIR_start=omp_get_wtime();
	FIR_para_INT_Mk3(input_data,spectra,coeff_mod,nTaps,nChannels,nBlocks,THREADS);
	FIR_end=omp_get_wtime();
	//********************************* FIR *********************************************
	
	if (check_errors) cout << "Error:" << FIR_check_uni(input_data,spectra,coeff_mod,nTaps,nChannels,nBlocks) << endl;
	
	//********************************* FFT *********************************************
	FFT_start=omp_get_wtime();
	FFT_mod(spectra,nTaps,nChannels,nBlocks);
	FFT_end=omp_get_wtime();
	//********************************* FFT *********************************************
	
	end=omp_get_wtime();
	//***************************** PERFORMANCE MEASUREMENTS *******************************
	*mtime_Decoupling=0;
	*mtime_FIR=FIR_end - FIR_start;
	*mtime_FFT=FFT_end - FFT_start;
	*mtime_Coupling=0;
	*mtime=end-start;
	
	delete[] input_data;
	delete[] coeff_mod;
	delete[] spectra;
}


void Test_intrinsic_Mk4(long int DATA_SIZE, int nTaps, int nChannels, double *mtime_Decoupling, double *mtime_FIR, double *mtime_FFT, double *mtime_Coupling, double *mtime){
	long int data_size=0,separated_data_size=0,f=0;
	double start, end, Decoupling_start, Decoupling_end, FIR_start, FIR_end, FFT_start, FFT_end, Coupling_start, Coupling_end;
	float *input_data;
	float *coeff_mod;
	float *spectra;
	int nBlocks;
	srand (time(NULL));
	
	data_size=DATA_SIZE*2;
	input_data = (float*)_mm_malloc( data_size*sizeof(float) ,64);
	spectra = (float*)_mm_malloc( data_size*sizeof(float) ,64);
	
	for(f=0;f<data_size;f++){
		input_data[f]=(float) (rand() % 10000 + 1)/1000.0;
	}
	
	separated_data_size=data_size/2;
	nBlocks=separated_data_size/nChannels; // determining number of blocks
	
	//***************************** PERFORMANCE MEASUREMENTS *******************************
	start=omp_get_wtime();
	
	//********************************* Window *********************************************
	coeff_mod = (float*)_mm_malloc( nChannels * nTaps * 2 * sizeof(float) ,64);
	Window_coefficients_mod(nTaps,nChannels,1,coeff_mod);
	//********************************* Window *********************************************
	
	//********************************* FIR *********************************************
	FIR_start=omp_get_wtime();
	FIR_para_INT_Mk4(input_data,spectra,coeff_mod,nTaps,nChannels,nBlocks,THREADS);
	FIR_end=omp_get_wtime();
	//********************************* FIR *********************************************
	
	if (check_errors) cout << "Error:" << FIR_check_uni(input_data,spectra,coeff_mod,nTaps,nChannels,nBlocks) << endl;
	
	//********************************* FFT *********************************************
	FFT_start=omp_get_wtime();
	FFT_mod(spectra,nTaps,nChannels,nBlocks);
	FFT_end=omp_get_wtime();
	//********************************* FFT *********************************************
	
	end=omp_get_wtime();
	//***************************** PERFORMANCE MEASUREMENTS *******************************
	*mtime_Decoupling=0;
	*mtime_FIR=FIR_end - FIR_start;
	*mtime_FFT=FFT_end - FFT_start;
	*mtime_Coupling=0;
	*mtime=end-start;
	
	_mm_free(input_data);
	_mm_free(coeff_mod);
	_mm_free(spectra);
}


void Test_intrinsic_Mk4_S(long int DATA_SIZE, int nTaps, int nChannels, double *mtime_Decoupling, double *mtime_FIR, double *mtime_FFT, double *mtime_Coupling, double *mtime){
	long int data_size=0,separated_data_size=0,f=0;
	double start, end, Decoupling_start, Decoupling_end, FIR_start, FIR_end, FFT_start, FFT_end, Coupling_start, Coupling_end;
	float *input_data;
	float *coeff_mod;
	float *spectra;
	int nBlocks;
	srand (time(NULL));
	
	data_size=DATA_SIZE*2;
	input_data = (float*)_mm_malloc( data_size*sizeof(float) ,64);
	spectra = (float*)_mm_malloc( data_size*sizeof(float) ,64);
	
	for(f=0;f<data_size;f++){
		input_data[f]=(float) (rand() % 10000 + 1)/1000.0;
	}
	
	separated_data_size=data_size/2;
	nBlocks=separated_data_size/nChannels; // determining number of blocks
	
	//***************************** PERFORMANCE MEASUREMENTS *******************************
	start=omp_get_wtime();
	
	//********************************* Window *********************************************
	coeff_mod = (float*)_mm_malloc( nChannels * nTaps * 2 * sizeof(float) ,64);
	Window_coefficients_mod(nTaps,nChannels,1,coeff_mod);
	//********************************* Window *********************************************
	
	//********************************* FIR *********************************************
	FIR_start=omp_get_wtime();
	FIR_serial_INT_Mk4(input_data,spectra,coeff_mod,nTaps,nChannels,nBlocks);
	FIR_end=omp_get_wtime();
	//********************************* FIR *********************************************
	
	if (check_errors) cout << "Error:" << FIR_check_uni(input_data,spectra,coeff_mod,nTaps,nChannels,nBlocks) << endl;
	
	//********************************* FFT *********************************************
	FFT_start=omp_get_wtime();
	FFT_mod(spectra,nTaps,nChannels,nBlocks);
	FFT_end=omp_get_wtime();
	//********************************* FFT *********************************************
	
	end=omp_get_wtime();
	//***************************** PERFORMANCE MEASUREMENTS *******************************
	*mtime_Decoupling=0;
	*mtime_FIR=FIR_end - FIR_start;
	*mtime_FFT=FFT_end - FFT_start;
	*mtime_Coupling=0;
	*mtime=end-start;
	
	_mm_free(input_data);
	_mm_free(coeff_mod);
	_mm_free(spectra);
}



void Test_intrinsic_Mk4_S_Mk1(long int DATA_SIZE, int nTaps, int nChannels, double *mtime_Decoupling, double *mtime_FIR, double *mtime_FFT, double *mtime_Coupling, double *mtime){
	long int data_size=0,separated_data_size=0,f=0;
	double start, end, Decoupling_start, Decoupling_end, FIR_start, FIR_end, FFT_start, FFT_end, Coupling_start, Coupling_end;
	float *input_data;
	float *coeff_mod;
	float *spectra;
	int nBlocks;
	srand (time(NULL));
	
	data_size=DATA_SIZE*2;
	input_data = (float*)_mm_malloc( data_size*sizeof(float) ,64);
	spectra = (float*)_mm_malloc( data_size*sizeof(float) ,64);
	
	for(f=0;f<data_size;f++){
		input_data[f]=(float) (rand() % 10000 + 1)/1000.0;
	}
	
	separated_data_size=data_size/2;
	nBlocks=separated_data_size/nChannels; // determining number of blocks
	
	//***************************** PERFORMANCE MEASUREMENTS *******************************
	start=omp_get_wtime();
	
	//********************************* Window *********************************************
	coeff_mod = (float*)_mm_malloc( nChannels * nTaps * 2 * sizeof(float) ,64);
	Window_coefficients_mod(nTaps,nChannels,1,coeff_mod);
	//********************************* Window *********************************************
	
	//********************************* FIR *********************************************
	FIR_start=omp_get_wtime();
	FIR_serial_Mk1_INT_Mk4(input_data,spectra,coeff_mod,nTaps,nChannels,nBlocks);
	FIR_end=omp_get_wtime();
	//********************************* FIR *********************************************
	
	if (check_errors) cout << "Error:" << FIR_check_uni(input_data,spectra,coeff_mod,nTaps,nChannels,nBlocks) << endl;
	
	//********************************* FFT *********************************************
	FFT_start=omp_get_wtime();
	FFT_mod(spectra,nTaps,nChannels,nBlocks);
	FFT_end=omp_get_wtime();
	//********************************* FFT *********************************************
	
	end=omp_get_wtime();
	//***************************** PERFORMANCE MEASUREMENTS *******************************
	*mtime_Decoupling=0;
	*mtime_FIR=FIR_end - FIR_start;
	*mtime_FFT=FFT_end - FFT_start;
	*mtime_Coupling=0;
	*mtime=end-start;
	
	_mm_free(input_data);
	_mm_free(coeff_mod);
	_mm_free(spectra);
}

void Test_intrinsic_Mk6(long int DATA_SIZE, int nTaps, int nChannels, double *mtime_Decoupling, double *mtime_FIR, double *mtime_FFT, double *mtime_Coupling, double *mtime){
	long int data_size=0,separated_data_size=0,f=0;
	double start, end, Decoupling_start, Decoupling_end, FIR_start, FIR_end, FFT_start, FFT_end, Coupling_start, Coupling_end;
	float *input_data;
	float *coeff_mod;
	float *spectra;
	int nBlocks;
	srand (time(NULL));
	
	data_size=DATA_SIZE*2;
	input_data = (float*)_mm_malloc( data_size*sizeof(float) ,64);
	spectra = (float*)_mm_malloc( data_size*sizeof(float) ,64);
	
	for(f=0;f<data_size;f++){
		input_data[f]=(float) (rand() % 10000 + 1)/1000.0;
	}
	
	separated_data_size=data_size/2;
	nBlocks=separated_data_size/nChannels; // determining number of blocks
	
	//***************************** PERFORMANCE MEASUREMENTS *******************************
	start=omp_get_wtime();
	
	//********************************* Window *********************************************
	coeff_mod = (float*)_mm_malloc( nChannels * nTaps * 2 * sizeof(float) ,64);
	Window_coefficients_mod(nTaps,nChannels,1,coeff_mod);
	//********************************* Window *********************************************
	
	//********************************* FIR *********************************************
	FIR_start=omp_get_wtime();
	FIR_para_INT_Mk6(input_data,spectra,coeff_mod,nTaps,nChannels,nBlocks,THREADS);
	FIR_end=omp_get_wtime();
	//********************************* FIR *********************************************
	
	if (check_errors) cout << "Error:" << FIR_check_uni(input_data,spectra,coeff_mod,nTaps,nChannels,nBlocks) << endl;
	
	//********************************* FFT *********************************************
	FFT_start=omp_get_wtime();
	FFT_mod(spectra,nTaps,nChannels,nBlocks);
	FFT_end=omp_get_wtime();
	//********************************* FFT *********************************************
	
	end=omp_get_wtime();
	//***************************** PERFORMANCE MEASUREMENTS *******************************
	*mtime_Decoupling=0;
	*mtime_FIR=FIR_end - FIR_start;
	*mtime_FFT=FFT_end - FFT_start;
	*mtime_Coupling=0;
	*mtime=end-start;
	
	_mm_free(input_data);
	_mm_free(coeff_mod);
	_mm_free(spectra);
}




int main(void) {
	THREADS=24;
	
	// -------> Time measuring
	double start, end, Decoupling_start, Decoupling_end, FIR_start, FIR_end, FFT_start, FFT_end, Coupling_start, Coupling_end;
    double mtime, mtime_Decoupling, mtime_FIR, mtime_FFT, mtime_Coupling;
	double sum_W, sum_Decoupling, sum_FIR, sum_FFT, sum_Coupling;
	double mean_W, mean_Decoupling, mean_FIR, mean_FFT, mean_Coupling;
    // -------> Time measuring

	// definition of constants
	int nTaps=8;
	int nChannels=512;
	const int nThreads=24;
	const int nRuns=10;
	
	long int DATA_BLOCKS=10000;
	long int DATA_SIZE=DATA_BLOCKS*nThreads*nChannels;
	long int MEMORY_READ=(2*DATA_SIZE + nTaps*nChannels*2*DATA_BLOCKS)*4;
	long int MEMORY_WRITE=(2*DATA_SIZE + 2*DATA_SIZE + nChannels)*4;
	long int FLOPS=nTaps*nChannels*DATA_BLOCKS*2;
	double flops,bandwidth;
	int count=0;
	
	bool debug=true;
	char str[200];
	
	ofstream PwI3_FILE,PwI4_FILE,SwI4_FILE,PwI6_FILE;
	PwI3_FILE.open("pwi3_C_perf.dat");PwI4_FILE.open("pwi4_C_perf.dat");
	SwI4_FILE.open("swi4_C_perf.dat");PwI6_FILE.open("pwi6_C_perf.dat");
	
	int data[5]={512,768,1024,2048,4096};

	for(int data_counter=0;data_counter<5;data_counter++){
		nChannels=data[data_counter];
		DATA_BLOCKS=2500;
		DATA_SIZE=DATA_BLOCKS*nThreads*nChannels;
		MEMORY_READ=2.0*DATA_SIZE*nTaps + 2.0*nTaps*DATA_SIZE;
		MEMORY_WRITE=2.0*DATA_SIZE;
		FLOPS=nTaps*DATA_SIZE*2.0 + nTaps*DATA_SIZE*2.0;
		cout << "-------------> Number of Channels: " << nChannels << endl;
		
		//********************************************************************************************
		//***************************** PERFORMANCE MEASUREMENTS *******************************
		//********************************************************************************************
		
		
		//**************
		// SERIAL Mk1 CODE WITH INTRINSIC Mk4
		//**************
		sum_Decoupling=0;sum_FIR=0;sum_FFT=0;sum_Coupling=0;sum_W=0;
		for(count=0;count<nRuns;count++){
			Test_intrinsic_Mk4_S_Mk1(DATA_SIZE,nTaps,nChannels,&mtime_Decoupling,&mtime_FIR,&mtime_FFT,&mtime_Coupling,&mtime);
			sum_Decoupling+=mtime_Decoupling;
			sum_FIR+=mtime_FIR;
			sum_FFT+=mtime_FFT;
			sum_Coupling+=mtime_Coupling;
			sum_W+=mtime;
			if(debug) cout << "SERIAL CODE Mk1 WITH INTRINSIC Mk4: FIR=" << mtime_FIR << " ;All=" << mtime << endl;
		}
		mean_W=(double) (sum_W/nRuns);mean_FIR=(double) (sum_FIR/nRuns);mean_FFT=(double) (sum_FFT/nRuns);
		
		bandwidth=((MEMORY_READ+MEMORY_WRITE)/mean_FIR)*4.0;
		flops=(FLOPS/mean_FIR);
		SwI4_FILE << data[data_counter] << " " << mean_FIR << " " << mean_W << " " << bandwidth << " " << flops << " " << mean_FFT << endl;
		
	
		//**************
		// PARALLEL CODE WITH INTRINSIC Mk3
		//**************
		sum_Decoupling=0;sum_FIR=0;sum_FFT=0;sum_Coupling=0;sum_W=0;
		for(count=0;count<nRuns;count++){
			Test_intrinsic_Mk3(DATA_SIZE,nTaps,nChannels,&mtime_Decoupling,&mtime_FIR,&mtime_FFT,&mtime_Coupling,&mtime);
			sum_Decoupling+=mtime_Decoupling;
			sum_FIR+=mtime_FIR;
			sum_FFT+=mtime_FFT;
			sum_Coupling+=mtime_Coupling;
			sum_W+=mtime;
			if(debug) cout << "PARALLEL CODE WITH INTRINSIC Mk3: FIR=" << mtime_FIR << " ;All=" << mtime << endl;
		}
		mean_W=(double) (sum_W/nRuns);mean_FIR=(double) (sum_FIR/nRuns);mean_FFT=(double) (sum_FFT/nRuns);
		
		bandwidth=((MEMORY_READ+MEMORY_WRITE)/mean_FIR)*4.0;
		flops=(FLOPS/mean_FIR);
		PwI3_FILE << data[data_counter] << " " << mean_FIR << " " << mean_W << " " << bandwidth << " " << flops << " " << mean_FFT << endl;
		
		
		/*
		//**************
		// SERIAL CODE WITH INTRINSIC Mk4
		//**************
		sum_Decoupling=0;sum_FIR=0;sum_FFT=0;sum_Coupling=0;sum_W=0;
		for(count=0;count<nRuns;count++){
			Test_intrinsic_Mk4_S(DATA_SIZE,nTaps,nChannels,&mtime_Decoupling,&mtime_FIR,&mtime_FFT,&mtime_Coupling,&mtime);
			sum_Decoupling+=mtime_Decoupling;
			sum_FIR+=mtime_FIR;
			sum_FFT+=mtime_FFT;
			sum_Coupling+=mtime_Coupling;
			sum_W+=mtime;
			if(debug) cout << "SERIAL CODE WITH INTRINSIC Mk4: FIR=" << mtime_FIR << " ;All=" << mtime << endl;
		}
		mean_W=(double) (sum_W/nRuns);mean_FIR=(double) (sum_FIR/nRuns);mean_FFT=(double) (sum_FFT/nRuns);
		
		bandwidth=((MEMORY_READ+MEMORY_WRITE)/mean_FIR)*4.0;
		flops=(FLOPS/mean_FIR);
		SwI4_FILE << data[data_counter] << " " << mean_FIR << " " << mean_W << " " << bandwidth << " " << flops << " " << mean_FFT << endl;
		*/
	
		
		//**************
		// PARALLEL CODE WITH INTRINSIC Mk4
		//**************
		sum_Decoupling=0;sum_FIR=0;sum_FFT=0;sum_Coupling=0;sum_W=0;
		for(count=0;count<nRuns;count++){
			Test_intrinsic_Mk4(DATA_SIZE,nTaps,nChannels,&mtime_Decoupling,&mtime_FIR,&mtime_FFT,&mtime_Coupling,&mtime);
			sum_Decoupling+=mtime_Decoupling;
			sum_FIR+=mtime_FIR;
			sum_FFT+=mtime_FFT;
			sum_Coupling+=mtime_Coupling;
			sum_W+=mtime;
			if(debug) cout << "PARALLEL CODE WITH INTRINSIC Mk4: FIR=" << mtime_FIR << " ;All=" << mtime << endl;
		}
		mean_W=(double) (sum_W/nRuns);mean_FIR=(double) (sum_FIR/nRuns);mean_FFT=(double) (sum_FFT/nRuns);
		
		bandwidth=((MEMORY_READ+MEMORY_WRITE)/mean_FIR)*4.0;
		flops=(FLOPS/mean_FIR);
		PwI4_FILE << data[data_counter] << " " << mean_FIR << " " << mean_W << " " << bandwidth << " " << flops << " " << mean_FFT << endl;		
		
		
		//**************
		// PARALLEL CODE WITH INTRINSIC Mk6
		//**************
		sum_Decoupling=0;sum_FIR=0;sum_FFT=0;sum_Coupling=0;sum_W=0;
		for(count=0;count<nRuns;count++){
			Test_intrinsic_Mk6(DATA_SIZE,nTaps,nChannels,&mtime_Decoupling,&mtime_FIR,&mtime_FFT,&mtime_Coupling,&mtime);
			sum_Decoupling+=mtime_Decoupling;
			sum_FIR+=mtime_FIR;
			sum_FFT+=mtime_FFT;
			sum_Coupling+=mtime_Coupling;
			sum_W+=mtime;
			if(debug) cout << "PARALLEL CODE WITH INTRINSIC Mk6: FIR=" << mtime_FIR << " ;All=" << mtime << endl;
		}
		mean_W=(double) (sum_W/nRuns);mean_FIR=(double) (sum_FIR/nRuns);mean_FFT=(double) (sum_FFT/nRuns);
		
		bandwidth=((MEMORY_READ+MEMORY_WRITE)/mean_FIR)*4.0;
		flops=(FLOPS/mean_FIR);
		PwI6_FILE << data[data_counter] << " " << mean_FIR << " " << mean_W << " " << bandwidth << " " << flops << " " << mean_FFT << endl;
		
		//********************************************************************************************
		//***************************** PERFORMANCE MEASUREMENTS *******************************
		//********************************************************************************************
		
	}
	PwI3_FILE.close();PwI4_FILE.close();SwI4_FILE.close();PwI6_FILE.close();
	
		
	cout << "Finished!" << endl;
	
    return EXIT_SUCCESS;
}