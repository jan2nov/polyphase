#include "timer.h"
//#include "utils.h"
#include "utils_cuda.h"
#include "data.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

typedef float2 Complex;

void reference_calculation(float2 *inputVals, float2 *outputVals, float *coeff, const int nChannels, unsigned int nBlocks, int nTaps);

void gpu_code(float2 *data, float2 *spectra, float *coeff, const int nChannels, unsigned int nBlocks, unsigned int filesize, int nThreads, int nTaps, int nStreams, int seg_blocks);

float reference_code(float2 *spectra_ref, float2 *spectra, int nChannels, unsigned int nTaps, unsigned int nBlocks);

int main(int argc, char **argv){
	
	int nTaps = 8;
	int nStreams = 2;
	int NUM_THREADS = 512;
	int seg_blocks = 5000;
	unsigned int data_size = 10000+nTaps-1;
	unsigned int nBlocks = 0;
	float error = 1.1f;
	bool debug=true;

	if (debug) printf("\t\tWelcome\n");

	Complex *h_signal, *h_spectra_pinned, *h_spectra_ref, *h_data_pinned;
	float *h_coeff;

	if (argc >= 2) NUM_THREADS = atof(argv[1]);
	if (argc >= 3) nStreams	  = (atof(argv[2]));
	if (argc >= 4) nTaps 	  = (atof(argv[3]));
	if (argc >= 5) seg_blocks = (atof(argv[4]));
	if (argc >= 6) data_size  = (atof(argv[5])+nTaps-1)*nChannels;

	nBlocks = (data_size+nTaps-1)/nChannels;nBlocks = data_size/nChannels;

	if (debug) printf("\nHost memory allocation...\t");
	checkCudaErrors(cudaMallocHost((void**)&h_spectra_pinned, data_size*sizeof(Complex)));
	checkCudaErrors(cudaMallocHost((void**)&h_data_pinned, data_size*sizeof(Complex)));
	h_signal 	= (Complex *)malloc(data_size*sizeof(Complex));
	h_spectra_ref = (Complex *)malloc(data_size*sizeof(Complex));
	h_coeff 	= (float *)malloc(nTaps*nChannels*sizeof(float));
	if (debug) printf("done.");

	if (debug) printf("\nHost memory memset...\t\t");
	memset(h_spectra_pinned, 0.0, sizeof(Complex)*data_size);	
	memset(h_spectra_ref, 0.0, sizeof(Complex)*data_size);	
	if (debug) printf("done.");

	if (debug) printf("\nLoad window coefficients...\t");
	//Load_window_data(h_coeff);
		for (int i = 0; i < nTaps*nChannels; i++)
			h_coeff[i] = rand() / (float)RAND_MAX;
	if (debug) printf("done.");


	if (debug) printf("\nRandom data set...\t\t");	
	srand(time(NULL));
	for (int i=0; i < (int)data_size; i++){
		h_signal[i].x = rand() / (float)RAND_MAX;
		h_signal[i].y = rand() / (float)RAND_MAX;
	}

	for (int i = 0; i < (int)data_size; i++){
		h_data_pinned[i] = h_signal[i];
	}
	if (debug) printf("done.");


	if (debug) printf("\nReference calculation...\t");
	reference_calculation(h_signal, h_spectra_ref, h_coeff, nChannels, nBlocks, nTaps);
	if (debug) printf("done.\n");
	
	//printf("CPU jedna %g druha %g", h_spectra_ref[3584], h_spectra_ref[7*512 + 259999].x);

	gpu_code(h_data_pinned, h_spectra_pinned, h_coeff, nChannels, nBlocks, data_size, NUM_THREADS, nTaps, nStreams, seg_blocks);	
	
	if (debug){
		error = reference_code(h_spectra_ref, h_spectra_pinned, nChannels, nTaps, nBlocks);
		printf( "error = %lf\n", error);
	}

	checkCudaErrors(cudaFreeHost(h_spectra_pinned));
	checkCudaErrors(cudaFreeHost(h_data_pinned));
	delete[] h_signal;
	delete[] h_spectra_ref;
	delete[] h_coeff;

}
