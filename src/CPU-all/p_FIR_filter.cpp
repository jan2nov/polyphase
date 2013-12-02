#include <omp.h>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <fftw3.h>
#include <string.h>
#include <xmmintrin.h>

#include <sys/time.h>
#include <unistd.h>

using namespace std;


long int File_size_byte(ifstream &FILEIN) {
	ifstream::pos_type size;
	FILEIN.seekg(0,ios::end);
	size=FILEIN.tellg();
	FILEIN.seekg(0,ios::beg);
	return((long int) size);
}

long int File_size_row(ifstream &FILEIN){
	std::size_t count=0;
	FILEIN.seekg(0,ios::beg);
	for(std::string line; std::getline(FILEIN, line); ++count){}
	return(count);
}


bool Load_data(ifstream &FILEIN, float *real, float *img, long int size){
	long int count=0;
	while (!FILEIN.eof()){
		FILEIN >> real[count];
		img[count]=0;
		count++;
	}
	return(1);
}


bool Save_data_spectra(char str[],float *spectra_real, float *spectra_img, unsigned int size, unsigned int nTaps, unsigned int nChannels){
	ofstream FILEOUT;
	FILEOUT.open(str);
	for (unsigned int c=0;c<size;c++){
		FILEOUT << spectra_real[c] << " " << spectra_img[c] << endl;
	}
	FILEOUT.close();
	return(1);
}

void Setup_buffer(float *w_buffer, float *real, float *img, unsigned int *oldesttap, int nChannels, int nTaps){
	*oldesttap=7;

	for(int i=0;i<nTaps;i++){
		for(int f=0;f<nChannels;f++){
			w_buffer[i*nChannels*2+f]=real[i*nChannels+f];
			w_buffer[i*nChannels*2+f+nChannels]=img[i*nChannels+f];
		}
	}

}

void Update_buffer(float *w_buffer, float *real, float *img, unsigned int *oldesttap, int nChannels, int nTaps){
	unsigned int itemp=*oldesttap;
	for(int f=0;f<nChannels;f++){
		w_buffer[itemp*nChannels*2+f]=real[f];
		w_buffer[itemp*nChannels*2+f+nChannels]=img[f];
	}
	itemp=(itemp+1)%nTaps;
	*oldesttap=itemp;
}

long int gettime(struct timeval start, struct timeval end){
	long int mtime,seconds,useconds;
    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
	return(mtime);
}

//****************************************************************************
//                           FIR filter with intrinsic
//****************************************************************************
void FIR_para_INT01(float *real, float *img, float *spectra_real, float *spectra_img, float *coeff, unsigned int nTaps, unsigned int nChannels, unsigned int nBlocks){
	unsigned int c,t,bl,tap,th_id,nthreads,block_step;
	unsigned int oldesttap;
	float *w_buffer;
	
	//registers
	unsigned int REG=128,FpR=4;
	__m128 i_real;
	__m128 i_img;
	__m128 i_coeff;
	__m128 i_spectra_real;
	__m128 i_spectra_img;
	__m128 i_temp;
	
	#pragma omp parallel shared(real,img,spectra_real,spectra_img) private(w_buffer,th_id,bl,t,c,oldesttap,tap,i_real,i_img,i_coeff,i_spectra_real,i_spectra_img,i_temp)
	{
		//omp_set_num_threads(48);
		w_buffer=new float[nChannels*nTaps*2];
		th_id = omp_get_thread_num();
		nthreads = omp_get_num_threads();
		block_step=nBlocks/nthreads;
		if(th_id==0) bl=0;
		else bl=(block_step*th_id)-(nTaps-1);
		
		Setup_buffer(w_buffer,&real[bl*nChannels],&img[bl*nChannels],&oldesttap,nChannels,nTaps);
		
		#pragma omp for schedule(static,block_step)
		for(bl=0;bl<nBlocks;bl++){
			
			Update_buffer(w_buffer,&real[bl*nChannels],&img[bl*nChannels],&oldesttap,nChannels,nTaps);
				
			for(c=0;c<REG;c++){
				// Zeroing!
				i_real=_mm_setzero_ps();
				i_img=_mm_setzero_ps();
				i_coeff=_mm_setzero_ps();
				i_spectra_real=_mm_setzero_ps();
				i_spectra_img=_mm_setzero_ps();
				i_temp=_mm_setzero_ps();
				
				for(t=0;t<nTaps;t++){
					tap=(oldesttap+t)%nTaps;
					//Loading data
					i_coeff=_mm_loadu_ps(&coeff[c*FpR+t*nChannels]);
					i_real=_mm_loadu_ps(&w_buffer[tap*2*nChannels+c*FpR]);
					i_img=_mm_loadu_ps(&w_buffer[tap*2*nChannels+c*FpR+nChannels]);

					i_temp=_mm_mul_ps(i_real,i_coeff);
					i_spectra_real=_mm_add_ps(i_temp,i_spectra_real);
					
					i_temp=_mm_mul_ps(i_img,i_coeff);
					i_spectra_img=_mm_add_ps(i_temp,i_spectra_img);
				
				} // for nTaps
				_mm_store_ps(&spectra_real[c*FpR+bl*nChannels],i_spectra_real);
				_mm_store_ps(&spectra_img[c*FpR+bl*nChannels],i_spectra_img);
			}// nChannels
		}// for nBlocks
		
		delete[] w_buffer;
		#pragma omp barrier
	}// parallel block
	
}


void FIR_para_FIR(float *real, float *img, float *spectra_real, float *spectra_img, float *coeff, unsigned int nTaps, unsigned int nChannels, unsigned int nBlocks){
	float *w_buffer;
	unsigned int c,t,bl,th_id,nthreads,block_step,tap,block;
	unsigned int oldesttap;
	
	#pragma omp parallel private(w_buffer,th_id,block,bl,t,c,oldesttap,tap) shared(real,img,spectra_real,spectra_img,coeff)
	{
		w_buffer=new float[nChannels*nTaps*2];
		th_id = omp_get_thread_num();
		nthreads = omp_get_num_threads();
		block_step=nBlocks/nthreads;
		if(th_id==0) block=0;
		else block=(block_step*th_id)-(nTaps-1);
		
		//******************** VALID ***********************
		Setup_buffer(w_buffer,&real[block*nChannels],&img[block*nChannels],&oldesttap,nChannels,nTaps);
		//******************** VALID ***********************
		#pragma omp for schedule(static,block_step)
		for(bl=0;bl<nBlocks;bl++){
					
			Update_buffer(w_buffer,&real[bl*nChannels],&img[bl*nChannels],&oldesttap,nChannels,nTaps);
			
			for(t=0;t<nTaps;t++){
				tap=(oldesttap+t)%nTaps;
				
				for(c=0;c<nChannels;c++){
					spectra_real[bl*nChannels+c]+=coeff[t*nChannels+c]*w_buffer[tap*2*nChannels+c];
					spectra_img[bl*nChannels+c]+=coeff[t*nChannels+c]*w_buffer[tap*2*nChannels+c+nChannels];
				} // nChannels
				
			} //nTaps
		}// nBlocks
		
		delete[] w_buffer;
		#pragma omp barrier
	}// parallel block
}

void FIR_serial(float *real, float *img, float *spectra_real, float *spectra_img, float *coeff, unsigned int nTaps, unsigned int nChannels, unsigned int nBlocks){
	float *w_buffer;
	unsigned int tap=0;
	unsigned int oldesttap;
	
	w_buffer=new float[nChannels*nTaps*2];

	Setup_buffer(w_buffer,real,img,&oldesttap,nChannels,nTaps);

	for(int bl=nTaps-1;bl<nBlocks;bl++){
		
		Update_buffer(w_buffer,&real[bl*nChannels],&img[bl*nChannels],&oldesttap,nChannels,nTaps);

		for(int t=0;t<nTaps;t++){
			tap=(oldesttap+t)%nTaps;
			for(int c=0;c<nChannels;c++){
				spectra_real[bl*nChannels+c]+=coeff[t*nChannels+c]*w_buffer[tap*2*nChannels+c];
				spectra_img[bl*nChannels+c]+=coeff[t*nChannels+c]*w_buffer[tap*2*nChannels+c+nChannels];
			}//nChannels
		}//nTaps
	} //nBlocks
}



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




int main(void) {
	// -------> Time measuring
	struct timeval start, end, FIR_start, FIR_end, FFT_start, FFT_end;
    long int mtime, mtime_FIR, mtime_FFT;
	double sum_W,sum_FIR,sum_FFT;
	double mean_W,mean_FIR,mean_FFT;
    // -------> Time measuring

	// definition of constants
	const int nTaps=8;
	const int nChannels=512;
	const int nThreads=24;
	
	long int DATA_BLOCKS=5000;
	long int DATA_SIZE=DATA_BLOCKS*nThreads*nChannels;
	long int MEMORY_READ=2*DATA_SIZE + nTaps*nChannels*2*DATA_BLOCKS;
	long int MEMORY_WRITE=2*DATA_SIZE + 2*DATA_SIZE + nChannels;
	long int FLOPS=nTaps*nChannels*DATA_BLOCKS*2;
	
	bool debug=true;
	char str[200];
	

	unsigned int nBlocks=0;
	long int filesize=0;
	long int old_filesize=0;

	float *real,*img;
	float *coeff;
	float *spectra_real,*spectra_img;
	bool flag=0;
	
	
	ofstream S_FILE,P_FILE,PwI_FILE;
	S_FILE.open("s_perf.dat");P_FILE.open("p_perf.dat");PwI_FILE.open("pwi_perf.dat");
	
	
	
	int data[5]={10,100,500,1000,5000};
		for(int data_counter=0;data_counter<5;data_counter++){
		DATA_BLOCKS=data[data_counter];
		DATA_SIZE=DATA_BLOCKS*nThreads*nChannels;
		MEMORY_READ=2*DATA_SIZE + nTaps*nChannels*2*DATA_BLOCKS;
		MEMORY_WRITE=2*DATA_SIZE + 2*DATA_SIZE + nChannels;
		FLOPS=nTaps*nChannels*DATA_BLOCKS*2;	
		cout << "Size of Data: " << DATA_BLOCKS << endl;

		//--------------> Loading data from file
		ifstream FILEIN;
		FILEIN.open("example1.dat",ios::in);
		if (!FILEIN.fail()){
			filesize=File_size_row(FILEIN);
			old_filesize=filesize;
			
			//******************** !! FOR TESTING !! ***********************
			old_filesize=filesize;
			filesize=DATA_SIZE;
			//******************** !! FOR TESTING !! ***********************
			
			real=new float[filesize];
			img=new float[filesize];
			
			FILEIN.clear();
			FILEIN.seekg(0,ios::beg);
			if (Load_data(FILEIN,real,img,old_filesize)) { // (- ! -)
				cout << "File opened and loaded" << endl;
			}
			
			//******************** !! FOR TESTING !! ***********************
			long int temp=0;
			if(filesize>old_filesize){
				for(long int f=0;f<(filesize-old_filesize);f++){
					temp=f%old_filesize;
					real[f+old_filesize]=real[temp];
					img[f+old_filesize]=img[temp];
				}
			}
			//******************** !! FOR TESTING !! ***********************
		}
		else {
			cout << "File not found" << endl;
			exit(2);
		}
		FILEIN.close();
		
		
		if(filesize%nChannels!=0){
			filesize=filesize-filesize%nChannels;
		}
		nBlocks=filesize/nChannels;
		

		spectra_real=new float[filesize];
		spectra_img=new float[filesize];
		memset(spectra_real,0,sizeof(float)*filesize);
		memset(spectra_img,0,sizeof(float)*filesize);
		
		
		//********************************************************************************************
		//***************************** PERFORMANCE MEASUREMENTS *******************************
		//********************************************************************************************
		int count=0;
		//********************************* Window *********************************************
		coeff = new float [nChannels * nTaps];
		Window_coefficients(nTaps,nChannels,1,coeff);
		//********************************* Window *********************************************
		
		//**************
		// SERIAL CODE
		//**************
		sum_W=0;sum_FIR=0;sum_FFT=0;
		for(count=0;count<100;count++){
			gettimeofday(&start, NULL);
			//********************************* FIR *********************************************
			gettimeofday(&FIR_start, NULL);
			FIR_serial(real,img,spectra_real,spectra_img,coeff,nTaps,nChannels,nBlocks);
			gettimeofday(&FIR_end, NULL);
			//********************************* FIR *********************************************
			
			//********************************* FFT *********************************************
			gettimeofday(&FFT_start, NULL);
			FFT(spectra_real,spectra_img,nTaps,nChannels,nBlocks);
			gettimeofday(&FFT_end, NULL);
			//********************************* FFT *********************************************
			gettimeofday(&end, NULL);
			
			mtime=gettime(start,end);
			mtime_FIR=gettime(FIR_start,FIR_end);
			mtime_FFT=gettime(FFT_start,FFT_end);
			sum_W+=mtime;sum_FIR+=mtime_FIR;sum_FFT+=mtime_FFT;
		}
		mean_W=(double) (sum_W/(count-1));mean_FIR=(double) (sum_FIR/(count-1));mean_FFT=(double) (sum_FFT/(count-1));
		cout << "Serial code" << endl;
		cout << "Mean times: W=" << mean_W << " sum_FIR=" << mean_FIR << " sum_FFT=" << mean_FFT << endl;
		cout << "Things: Bandwidth=" << ((MEMORY_READ+MEMORY_WRITE)/mean_FIR)*1000.0 << " Flops=" << (FLOPS/mean_FIR)*1000.0 << endl;
		S_FILE << data[data_counter] << " " << mean_FIR << " " << ((MEMORY_READ+MEMORY_WRITE)/mean_FIR)*1000.0 << " " << (FLOPS/mean_FIR)*1000.0 << endl;
		
		
		//**************
		// PARALLEL CODE
		//**************
		sum_W=0;sum_FIR=0;sum_FFT=0;
		for(count=0;count<100;count++){
			gettimeofday(&start, NULL);
			//********************************* FIR *********************************************
			gettimeofday(&FIR_start, NULL);
			FIR_para_FIR(real,img,spectra_real,spectra_img,coeff,nTaps,nChannels,nBlocks);
			gettimeofday(&FIR_end, NULL);
			//********************************* FIR *********************************************
			
			//********************************* FFT *********************************************
			gettimeofday(&FFT_start, NULL);
			FFT(spectra_real,spectra_img,nTaps,nChannels,nBlocks);
			gettimeofday(&FFT_end, NULL);
			//********************************* FFT *********************************************
			gettimeofday(&end, NULL);
			
			mtime=gettime(start,end);
			mtime_FIR=gettime(FIR_start,FIR_end);
			mtime_FFT=gettime(FFT_start,FFT_end);
			sum_W+=mtime;sum_FIR+=mtime_FIR;sum_FFT+=mtime_FFT;
		}
		mean_W=(double) (sum_W/(count-1));mean_FIR=(double) (sum_FIR/(count-1));mean_FFT=(double) (sum_FFT/(count-1));
		cout << "Parallel code" << endl;
		cout << "Mean times: W=" << mean_W << " sum_FIR=" << mean_FIR << " sum_FFT=" << mean_FFT << endl;
		cout << "Things: Bandwidth=" << ((MEMORY_READ+MEMORY_WRITE)/mean_FIR)*1000.0 << " Flops=" << (FLOPS/mean_FIR)*1000.0 << endl;
		P_FILE << data[data_counter] << " " << mean_FIR << " " << ((MEMORY_READ+MEMORY_WRITE)/mean_FIR)*1000.0 << " " << (FLOPS/mean_FIR)*1000.0 << endl;
		
		
		//**************
		// PARALLEL CODE WITH INTRINSIC
		//**************
		sum_W=0;sum_FIR=0;sum_FFT=0;
		for(count=0;count<100;count++){
			gettimeofday(&start, NULL);
			//********************************* FIR *********************************************
			gettimeofday(&FIR_start, NULL);
			FIR_para_INT01(real,img,spectra_real,spectra_img,coeff,nTaps,nChannels,nBlocks);
			gettimeofday(&FIR_end, NULL);
			//********************************* FIR *********************************************
			
			//********************************* FFT *********************************************
			gettimeofday(&FFT_start, NULL);
			FFT(spectra_real,spectra_img,nTaps,nChannels,nBlocks);
			gettimeofday(&FFT_end, NULL);
			//********************************* FFT *********************************************
			gettimeofday(&end, NULL);
			
			mtime=gettime(start,end);
			mtime_FIR=gettime(FIR_start,FIR_end);
			mtime_FFT=gettime(FFT_start,FFT_end);
			sum_W+=mtime;sum_FIR+=mtime_FIR;sum_FFT+=mtime_FFT;
		}
		mean_W=(double) (sum_W/(count-1));mean_FIR=(double) (sum_FIR/(count-1));mean_FFT=(double) (sum_FFT/(count-1));
		cout << "Parallel code with intrinsic" << endl;
		cout << "Mean times: W=" << mean_W << " sum_FIR=" << mean_FIR << " sum_FFT=" << mean_FFT << endl;
		cout << "Things: Bandwidth=" << ((MEMORY_READ+MEMORY_WRITE)/mean_FIR)*1000.0 << " Flops=" << (FLOPS/mean_FIR)*1000.0 << endl;
		PwI_FILE << data[data_counter] << " " << mean_FIR << " " << ((MEMORY_READ+MEMORY_WRITE)/mean_FIR)*1000.0 << " " << (FLOPS/mean_FIR)*1000.0 << endl;
		
		
		//********************************************************************************************
		//***************************** PERFORMANCE MEASUREMENTS *******************************
		//********************************************************************************************
		
		
		delete[] coeff;
		delete[] real;
		delete[] img;
		delete[] spectra_real;
		delete[] spectra_img;
	}
	S_FILE.close();P_FILE.close();PwI_FILE.close();
	
	
	
	
	/*
	// -------> Time measuring
	mtime=gettime(start,end);
	mtime_FIR=gettime(FIR_start,FIR_end);
	mtime_FFT=gettime(FFT_start,FFT_end);
	printf("FIR filter time: %ld milliseconds\n", mtime_FIR);
	printf("FFT time: %ld milliseconds\n", mtime_FFT);
	printf("Elapsed time: %ld milliseconds\n", mtime);
	// -------> Time measuring
	*/
    return EXIT_SUCCESS;
}