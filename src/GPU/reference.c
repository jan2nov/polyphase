#include <stdio.h>
#include "data.h"
//#include <cuda.h>
#include <cufft.h>
#include <stdlib.h>
#include <fftw3.h>

float reference_code(float2 *spectra_ref, 
		     float2 *spectra, 
		     const int nChannels, 
		     unsigned int nTaps, 
		     unsigned int nBlocks){

	double diff = 0.0f, error_norm = 0.0f;
	for (int j = 0; j < (int)nBlocks - (int)nTaps + 1; j++){
	  for (int i = 0; i < nChannels; i++) {
	    diff = spectra_ref[i + (nTaps-1)*nChannels + j*nChannels].x - spectra[j*nChannels + i].x;
	      error_norm += diff * diff;
	      //if (diff != 0.0) printf("%i %g %g\n", i +j*nChannels  , spectra_ref[i + 7*nChannels+ j*nChannels ].x, spectra[j*nChannels + i].x);
	      diff = spectra_ref[i + (nTaps-1)*nChannels + j*nChannels].y - spectra[j*nChannels + i].y;
	      error_norm += diff * diff;
	  }
	}
	error_norm = (float)sqrt((double)error_norm);

	return error_norm;
}

void Setup_buffer(float *buffer, 
		  float2 *data, 
		  unsigned int *oldesttap, 
		  int nChannels, 
		  int nTaps){

	*oldesttap = nTaps - 1;

	for(int i = 0; i < nTaps; i++){
		for(int f = 0; f < nChannels; f++){
			// structure in memory [nchannel_real, nchannel_img] per one tap
			buffer[i*nChannels*2+f]           = data[i*nChannels+f].x;
			buffer[i*nChannels*2+f+nChannels] = data[i*nChannels+f].y;
		}
	}
}

void Update_buffer(float *w_buffer, 
		   float2 *data, 
		   unsigned int *oldesttap, 
		   int nChannels, 
		   int nTaps){

	unsigned int itemp=*oldesttap;

	for(int f=0; f<nChannels; f++){
		w_buffer[itemp*nChannels*2+f]		= (float)data[f].x;
		w_buffer[itemp*nChannels*2+f+nChannels] = data[f].y;
	}

	itemp=(itemp+1)%nTaps;
	*oldesttap=itemp;
}

void Fir_cpu(float *w_buffer, 
	     float2 *data, 
	     unsigned int *oldesttap, 
	     int nChannels, 
	     int nTaps, 
	     int nBlocks, 
         float *coeff, 
 	     float2 *spectra){

	unsigned int tap = 0;

	for(int bl=nTaps - 1; bl < nBlocks; bl++){
		Update_buffer(w_buffer, &data[bl*nChannels], oldesttap, nChannels, nTaps);
			for(int t=0; t<nTaps; t++){
				tap=(*oldesttap+t)%(nTaps);	
				//printf("%d\n", tap);
				for(int c=0; c < nChannels; c++){ 
				  spectra[bl*nChannels+c].x += coeff[t*nChannels+c]*w_buffer[tap*2*nChannels+c];
				  spectra[bl*nChannels+c].y += coeff[t*nChannels+c]*w_buffer[tap*2*nChannels+c+nChannels];
				}
			}
	}
}

void  reference_calculation(float2 *inputVals, 
			    float2 *outputVals,
			    float *coeff, 
			    int nChannels, 
			    unsigned int nBlocks,
			    int nTaps){

	float *w_buffer;
	unsigned int oldesttap;

	w_buffer = (float *)malloc(sizeof(float)*nChannels*nTaps*2);

	Setup_buffer(w_buffer, inputVals, &oldesttap, nChannels, nTaps);
	Fir_cpu(w_buffer, inputVals, &oldesttap, nChannels, nTaps, nBlocks, coeff, outputVals);

	/*	fftwf_plan p;
		fftwf_complex *in, *out;

		in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChannels);
		out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChannels);
		p = fftwf_plan_dft_1d(nChannels, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		if (!p){
			exit(101);
		}

		for(int bl=nTaps-1;bl<(int)nBlocks;bl++){
			// I need better memory arrangement
			for(int f=0; f<nChannels; f++){
				in[f][0] = outputVals[bl*nChannels+f].x;
				in[f][1] = outputVals[bl*nChannels+f].y;
			}
			fftwf_execute(p);
			for(int f=0;f<nChannels;f++){
				outputVals[bl*nChannels+f].x = out[f][0];
				outputVals[bl*nChannels+f].y = out[f][1];
			}
		}

	fftwf_destroy_plan(p);
	fftwf_free(in);
	fftwf_free(out);*/
	delete[] w_buffer;
}
