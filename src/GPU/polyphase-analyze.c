#include "timer.h"
#include "utils_cuda.h"
#include "data.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

typedef float2 Complex;

void reference_calculation(float2 *inputVals, float2 *outputVals, float *coeff, const int nChannels, unsigned int nBlocks, int nTaps);

void gpu_code(float2 *data_in, float2 *spectra, float *coeff, const int nChannels, unsigned int nBlocks, unsigned int filesize, const int nThreads, int blocks_y, const int nTaps);

float reference_code(float2 *spectra_ref, float2 *spectra, int nChannels, unsigned int nTaps, unsigned int nBlocks);

int main(int argc, char **argv){
	
	int nTaps = 8;
	int num_threads = 512;
	int NUM_BLOCKS = 1;
	unsigned int data_size = 10000+nTaps-1;
	unsigned int nBlocks = 0;
	float error = 1.1f;
	bool debug=true;

	if (debug) printf("\t\tWelcome\n");

	Complex *h_signal, *h_spectra, *h_spectra_ref;
	float *h_coeff;

	if (argc >= 2) NUM_BLOCKS  = atof(argv[1]);
	if (argc >= 3) num_threads = atof(argv[2]);
	if (argc >= 4) nTaps 	   = (atof(argv[3]));
	if (argc >= 5) data_size   = (atof(argv[4])+nTaps-1)*nChannels;

	nBlocks = (data_size)/nChannels;

	if (debug) printf("\nHost memory allocation...\t");
		h_signal 	= (Complex *)malloc(data_size*sizeof(Complex));
		h_spectra 	= (Complex *)malloc((data_size)*sizeof(Complex));
		h_spectra_ref = (Complex *)malloc((data_size)*sizeof(Complex));
		h_coeff 	= (float *)malloc(nTaps*nChannels*sizeof(float));
	if (debug) printf("done.");

	if (debug) printf("\nHost memory memset...\t\t");
		memset(h_spectra, 0.0, sizeof(Complex)*(data_size));	
		memset(h_spectra_ref, 0.0, sizeof(Complex)*(data_size));	
	if (debug) printf("done.");

	if (debug) printf("\nLoad window coefficients...\t");
	//Load_window_data(h_coeff);
		srand(time(NULL));
		for (int i = 0; i < nTaps*nChannels; i++)		
			//h_coeff[i] = 1.0;
			h_coeff[i] = rand() / (float)RAND_MAX;
	if (debug) printf("done.");


	if (debug) printf("\nRandom data set...\t\t");	
		srand(time(NULL));
		for (int i=0; i < (int)data_size; i++){
			h_signal[i].x = rand() / (float)RAND_MAX;
			h_signal[i].y = rand() / (float)RAND_MAX;
		}

	if (debug) printf("done.");

	//printf("\n%g %g\n", h_real[0], h_real[65536*512]);

	if (debug) printf("\nReference calculation...\t");
		reference_calculation(h_signal, h_spectra_ref, h_coeff, nChannels, nBlocks, nTaps);
	if (debug) printf("done.\n");

	gpu_code(h_signal, h_spectra, h_coeff, nChannels, nBlocks, data_size, num_threads, NUM_BLOCKS, nTaps);	
	
	if (debug){
		error = reference_code(h_spectra_ref, h_spectra, nChannels, nTaps, nBlocks);
		printf( "error = %lf\n", error);
	}

	delete[] h_signal;
	delete[] h_spectra;
	delete[] h_spectra_ref;
	delete[] h_coeff;
	
	checkCudaErrors(cudaDeviceReset());
	
	printf("\nThat's All folks!\n");

	return 0;
}
