#include "timer.h"
//#include "utils.h"
#include "utils_file.h"
#include "data.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

typedef float2 Complex;

void reference_calculation(float2 *inputVals, float2 *outputVals, float *coeff, const int nChannels, unsigned int nBlocks);

void gpu_code(float *real, float *img, float2 *spectra, float *coeff, const int nChannels, unsigned int nBlocks, unsigned int filesize, int blocks_y);

float reference_code(float2 *spectra_ref, float2 *spectra, int nChannels, unsigned int nBlocks);

int main(int argc, char **argv){
	
	int NUM_BLOCKS = 1;
	unsigned int data_size = 512512;
	unsigned int nBlocks = 0;
	float error = 1.1f;
	bool debug=false;

	Complex *h_signal, *h_spectra, *h_spectra_ref;
	float *h_coeff, *h_real, *h_img;

	if (argc >= 2) NUM_BLOCKS = atof(argv[1]);
	if (argc >= 3) data_size  = atof(argv[2])*nChannels;

	nBlocks = data_size/nChannels;

	h_signal 	= (Complex *)malloc(data_size*sizeof(Complex));
	h_real	 	= (float *)malloc(data_size*sizeof(float));
	h_img	 	= (float *)malloc(data_size*sizeof(float));
	h_spectra 	= (Complex *)malloc(data_size*sizeof(Complex));
	h_spectra_ref 	= (Complex *)malloc(data_size*sizeof(Complex));
	h_coeff 	= (float *)malloc(nTaps*nChannels*sizeof(float));
	
	memset(h_spectra, 0.0, sizeof(Complex)*data_size);	
	memset(h_spectra_ref, 0.0, sizeof(Complex)*data_size);	

	Load_window_data(h_coeff);

	srand(time(NULL));
	for (int i=0; i < (int)data_size; i++){
		h_signal[i].x = rand() / (float)RAND_MAX;
		h_signal[i].y = rand() / (float)RAND_MAX;
	}

	for (int i = 0; i < (int)data_size; i++){
		h_real[i] = h_signal[i].x;
		h_img[i]  = h_signal[i].y;
	}

	reference_calculation(h_signal, h_spectra_ref, h_coeff, nChannels, nBlocks);

	gpu_code(h_real, h_img, h_spectra, h_coeff, nChannels, nBlocks, data_size, NUM_BLOCKS);

	if (debug){
		error = reference_code(h_spectra_ref, h_spectra, nChannels, nTaps);
		printf( "error = %g\n", error);
	}

	delete[] h_signal;
	delete[] h_real;
	delete[] h_img;
	delete[] h_spectra;
	delete[] h_spectra_ref;
	delete[] h_coeff;

}
