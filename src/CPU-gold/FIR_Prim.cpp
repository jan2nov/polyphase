//TODO! if I will keep the current buffer arangement with blocks of size nTaps I need to change coeff to be more compatible with the chosen scheme
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
#include <math.h>
#include <fftw3.h>

#include "Filter_window.h"

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
	for (unsigned int c=(nTaps-1)*nChannels;c<size;c++){
		FILEOUT << spectra_real[c] << " " << spectra_img[c] << endl;
	}
	FILEOUT.close();
	return(1);
}

bool Save_window_data(float *coeff, unsigned int nChannels, unsigned int nTaps){
	ofstream FILEOUT;
	FILEOUT.open("window.tap");

	for(unsigned int t=0;t<nTaps;t++){
		for (unsigned int c=0;c<nChannels;c++){
			FILEOUT << coeff[t*nChannels+c] << endl;
		}
	}

	FILEOUT.close();
	return (true);
}

bool Load_window_data(ifstream &FILEIN, float *coeff){
	unsigned int count=0;
	ifstream FILEIN;
	FILEIN.open("window.tap",ios::in);
	while (!FILEIN.eof()){
		FILEIN >> coeff[count];
		count++;
	}
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

int main(void) {
	bool debug=true;
	char str[200];
	
	// definition of constants
	const int nTaps=8;
	const int nChannels=512;
	//const int register_size=256;//CPU
	//const int in_buffer_size=register_size/32;	
	//const int w_buffer_size=nChannels*nTaps*3;

	unsigned int nBlocks=0;
	unsigned int oldesttap;
	long int filesize=0;

	float *real,*img;
	float *coeff;
	float *spectra_real,*spectra_img;
	float *w_buffer;
	bool flag=0;

	w_buffer=new float[nChannels*nTaps*2];
	
	//--------------> Loading data from file
	ifstream FILEIN;
	FILEIN.open("example3.dat",ios::in);
	if (!FILEIN.fail()){
		filesize=File_size_row(FILEIN);
		
		real=new float[filesize];
		img=new float[filesize];
		spectra_real=new float[filesize];
		spectra_img=new float[filesize];
		for(long int f=0;f<filesize;f++){
			spectra_real[f]=0;
			spectra_img[f]=0;
		}

		FILEIN.clear();
		FILEIN.seekg(0,ios::beg);
		flag=Load_data(FILEIN,real,img,filesize);
		if (flag) { 
			cout << "File opened and loaded" << endl;
		}
	}
	else {
		cout << "File not found" << endl;
		exit(2);
	}
	FILEIN.close();
	
	if(filesize%nChannels!=0){
		filesize=filesize-filesize%nChannels;
		if(debug) cout << filesize << " - " << filesize/512 << endl;
	}
	nBlocks=filesize/nChannels;
	if(debug) cout << "Number of rows: " << filesize << endl;
	//-----------------------------------------------------------------------

	//------------------> Getting window coefficients
	coeff = new float [nChannels * nTaps];
	//Creating filter window
	Window_coefficients(nTaps,nChannels,1,coeff);
	Save_window_data(coeff,nChannels,nTaps);
	//-----------------------------------------------------------------------
	
	
	//****************************************************************************
	//                           Primitive FIR filter
	//****************************************************************************
	double temp1=0,temp2=0,temp3=0;
	unsigned int tap=0;

	Setup_buffer(w_buffer,real,img,&oldesttap,nChannels,nTaps);

	for(int bl=nTaps-1;bl<nBlocks;bl++){
		Update_buffer(w_buffer,&real[bl*nChannels],&img[bl*nChannels],&oldesttap,nChannels,nTaps);
		for(int t=0;t<nTaps;t++){
			tap=(oldesttap+t)%(nTaps-1);
			for(int c=0;c<nChannels;c++){
				temp1=coeff[t*nChannels+c];
				temp2=w_buffer[tap*2*nChannels+c];
				spectra_real[bl*nChannels+c]+=coeff[t*nChannels+c]*w_buffer[tap*2*nChannels+c];
				spectra_img[bl*nChannels+c]+=coeff[t*nChannels+c]*w_buffer[tap*2*nChannels+c+nChannels];
			}
		}
	}

	sprintf(str,"outtest_wth_DFT.dat");
	Save_data_spectra(str,spectra_real,spectra_img,filesize,nTaps,nChannels);


	//****************************************************************************
	//                           Primitive FFT using FFTW
	//****************************************************************************
	int N=nChannels;
	fftwf_complex *in, *out;
	fftwf_plan p;
	
	in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
	out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
	p = fftwf_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	if (!p){
		exit(101);
	}
	for(int bl=nTaps-1;bl<nBlocks;bl++){
		// I need better memory arrangement
		for(int f=0;f<N;f++){
			in[f][0]=spectra_real[bl*nChannels+f];
			in[f][1]=spectra_img[bl*nChannels+f];
		}
		fftwf_execute(p);
		for(int f=0;f<N;f++){
			spectra_real[bl*nChannels+f]=out[f][0];
			spectra_img[bl*nChannels+f]=out[f][1];
		}
	}
	fftwf_destroy_plan(p);
	fftwf_free(in);
	fftwf_free(out);

	sprintf(str,"outtest_w_DFT.dat");
	Save_data_spectra(str,spectra_real,spectra_img,filesize,nTaps,nChannels);






	delete[] coeff;
	delete[] real;
	delete[] img;
	delete[] spectra_real;
	delete[] spectra_img;
	delete[] w_buffer;

    #pragma omp parallel
    printf("Hello, world.\n");
    system("PAUSE");
    return EXIT_SUCCESS;
}

