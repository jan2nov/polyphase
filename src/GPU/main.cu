#include <stdio.h>
#include <cuda.h>
#include <cufft.h>
#include <fftw3.h>
#include <fstream>
#include <iostream>
#include <cutil_inline.h>
#include <time.h>

using namespace std;

#include "utils_file.h"
#include "utils_reference.h"
#include "utils_cuda.h"
#include "timer.h"

// complex data type
typedef float2 Complex;
__constant__ float C_coeff[4096];

void Setup_buffer(float *w_buffer, Complex *data, unsigned int *oldesttap, int nChannels, int nTaps){
	*oldesttap = 7;

	for(int i = 0; i < nTaps; i++){
		for(int f = 0; f < nChannels; f++){
			// structure in memory [nchannel_real, nchannel_img] per one tap
			w_buffer[i*nChannels*2+f]           = data[i*nChannels+f].x;
			w_buffer[i*nChannels*2+f+nChannels] = data[i*nChannels+f].y;
		}
	}
}

void Update_buffer(float *w_buffer, Complex *data, unsigned int *oldesttap, int nChannels, int nTaps){
	unsigned int itemp=*oldesttap;
	for(int f=0; f<nChannels; f++){
		w_buffer[itemp*nChannels*2+f]=data[f].x;
		w_buffer[itemp*nChannels*2+f+nChannels]=data[f].y;
	}
	itemp=(itemp+1)%nTaps;
	*oldesttap=itemp;
}

void Fir_cpu(float *w_buffer, Complex *data, unsigned int *oldesttap, int nChannels, int nTaps, int nBlocks, float *coeff, Complex *spectra){

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

__global__ void Fir_old(float *d_real,float *d_img, float *d_coeff, int nTaps, int nChannels, int nBlocks, float *d_spectra_real, float *d_spectra_img) {
	int oldesttap = 0;
	unsigned int tap = 0;
	int c = threadIdx.x;

	for(int bl = nTaps-1; bl < nBlocks; bl++){
		for(int t=0;t<nTaps;t++){
			tap=(oldesttap+t)%(nTaps);
			d_spectra_real[c] += d_coeff[t*nChannels+c]*d_real[tap*nChannels+c];
			d_spectra_img[c]  += d_coeff[t*nChannels+c]*d_img[tap*nChannels+c];
		}
	}
}

__global__ void Fir(float *d_signal_real, float *d_signal_img, int nTaps, int nChannels, Complex *spectra)
{
	int tx = threadIdx.x + blockDim.x*blockIdx.y;
	int index = nChannels*blockIdx.x + tx;
	int i, i_coeff, i_data;
	float2 local_spectra = {0.0f,0.0f};
	//printf("%g\n", C_coeff[0]);

	for(int t=0;t<nTaps;t++){
	  i = t*nChannels;
	  i_coeff = i + tx;
	  i_data = i + index;
	  local_spectra.x += C_coeff[i_coeff]*d_signal_real[i_data];
	  // spectra[index].x += C_coeff[t*nChannels + tx];//*d_signal_real[t*nChannels + index];
	  //	spectra[index].y += C_coeff[t*nChannels + tx]*d_signal_img[t*nChannels + index];
	  local_spectra.y += C_coeff[i_coeff]*d_signal_img[i_data];
	}
		spectra[index]=local_spectra;
	
}

__global__ void Fir_01(float *d_signal_real, float *d_signal_img, float *coeff, int nTaps, int nChannels, Complex *spectra)
{
	extern __shared__ float in_data[];

	int tx = threadIdx.x + blockDim.x*blockIdx.y;
	int index = blockIdx.x*nChannels + tx;

	float *real = &in_data[0];
	float *img  = &in_data[nTaps*nChannels];
	float sum_real = 0, sum_img = 0;

	//#pragma unroll 10
	for(int t=0; t < nTaps; t++){
		//real[index] = d_signal_real[index];
		real[tx + t*nChannels] = d_signal_real[index + t*nChannels];
		img[tx + t*nChannels] = d_signal_img[t*nChannels + index];
	}
	__syncthreads();

	for(int t=0;t<nTaps;t++){
		sum_real += C_coeff[t*nChannels + tx]*real[t*nChannels + tx];
		sum_img  += C_coeff[t*nChannels + tx]*img[t*nChannels + tx];
	}
//	if (index == 0) printf("Cislo: %g %g %d\n",in_data[4096], d_signal_img[0], blockDim.x);

	spectra[index].x = sum_real;
	spectra[index].y = sum_img;
	return;
}

__global__ void test_01(float *in, float*out)
{
	int tx = threadIdx.x;
//	int bx = blockIdx.x;
//	printf("%g",in[32]);
	for (int i = 0; i < gridDim.x;i++){
		out[tx] += in[i*blockDim.x + tx];
		__syncthreads();
	}

	//out[tx] = blockDim.x;
}

__global__ void test(Complex *data){
	__shared__ float sdata[1024];

	float* ptr = (float *) &data[0];
	float2 *o_data = (float2 *) sdata;

	int index = blockIdx.x*blockDim.x + threadIdx.x;

	sdata[index] = ptr[index];
	sdata[index + blockDim.x] = ptr[blockDim.x + index];
	__syncthreads();

	o_data[index].x += 2;
	__syncthreads();

	ptr[index] = sdata[index];
	ptr[index + blockDim.x] = sdata[index + blockIdx.x];
}

__global__ void complex2real(Complex *data, float *data_real, float *data_img){

	int index = blockIdx.x*blockDim.x + threadIdx.x;
	float2 store = {0.0f,0.0f};

	store = data[index];

	data_real[index] = store.x;
	data_img[index]  = store.y;

	return;
}

int main(int argc,char **argv)
{
		//for debug use only
		bool debug=true;
		//set some variables
		srand(time(NULL));
		GpuTimer timer;
		int NUM_BLOCKS = 1;
		int BLOCK_WIDTH = 512;
		int TEST = 1; 
		unsigned int nBlocks = 0;
		Complex *h_signal, *h_spectra;
		const int nTaps = 8;
		const int nChannels = 512;

		int N = 512512;

		//      		cudaSetDevice(1);

		if (argc >= 2) NUM_BLOCKS = atof(argv[1]);
		if (argc >= 3) N = atof(argv[2])*nChannels;
		if (argc == 4) TEST = atof(argv[3]);
	        printf("N: %i\n", N);
		//float *d_data, *d_odata;
		h_signal = new Complex[N];
		h_spectra= new Complex[N];
		//		for (int j = 0; j<M; j++){
			for (int i=0; i < N; i++){
			  h_signal[i].x = rand() / (float)RAND_MAX;
			  h_signal[i].y = 0.0f;
			}
			//}
			//		checkCudaErrors(cudaMalloc((void **) &d_data,  N*sizeof(float)));
			//		checkCudaErrors(cudaMalloc((void **) &d_odata,  N*sizeof(Complex)));
			//		checkCudaErrors(cudaMemset(d_odata,  0.0, N*sizeof(Complex)));
			//checkCudaErrors(cudaMemcpy(d_data, h_data, N*sizeof(float), cudaMemcpyHostToDevice));

      		if (debug) printf("Num of thread blocks: %d\n", NUM_BLOCKS);
		if (debug) printf("x: %g y: %g\n", h_signal[0].x, h_signal[0].y);
		//-----------------------------------------------------------------------
		//---------------------------> load data
		ifstream FILEIN;
		long int filesize = 0;
		bool flag = 0;
		/*
		FILEIN.open("fakedata1.dat",ios::in);
		//FILEIN.open("example1.tap",ios::in);
			if (!FILEIN.fail()){
				filesize = File_size_row(FILEIN);

				h_signal  = new Complex[filesize];
				h_spectra = new Complex[filesize];
				// cleaning the variables
				for(long int f = 0; f < filesize; f++){
					h_spectra[f].x = 0;
					h_spectra[f].y = 0;
				}

				FILEIN.clear();
				FILEIN.seekg(0,ios::beg);

				//load data to h_real and h_img variables
				flag = Load_data(FILEIN, h_signal, filesize);
				if (flag){
				  if (debug) cout << "File opened and loaded" << endl;
				}
			}
			else {
				cout << "File not found" << endl;
				exit(2);
			}
			FILEIN.close();
		*/
		filesize = N;
		if(filesize%nChannels != 0){
				filesize = filesize - filesize%nChannels;
				if(debug) cout << filesize << " - " << filesize / 512 << endl;
		}

		nBlocks = filesize/nChannels;
		if(debug) cout << "Number of rows: " << filesize << endl;
		if(debug) cout << "Number of Blocks: " << nBlocks  << endl;
		if(debug) cout << " " << endl;

		//-----------------------------------------------------------------------
	    //------------------> Getting window coefficients
		float *h_coeff;
		h_coeff = new float [nChannels * nTaps];

		//Creating filter window
		Load_window_data(h_coeff);


/////////////////////////////////////////////////////////////////////////////
// reference checking
/////////////////////////////////////////////////////////////////////////////

				//Allocate data for reference
				Complex *h_spectra_ref;
				unsigned int oldesttap;
				float *w_buffer;
				float error;

				h_spectra_ref = new Complex[filesize];
				w_buffer      = new float[nChannels*nTaps*2];

				memset(h_spectra_ref, 0.0, sizeof(Complex)*filesize);
				memset(h_spectra, 0.0, sizeof(Complex)*filesize);
				// compute CPU FIR code
				Setup_buffer(w_buffer, h_signal, &oldesttap, nChannels, nTaps);
				Fir_cpu(w_buffer, h_signal, &oldesttap, nChannels, nTaps, nBlocks, h_coeff, h_spectra_ref);
				//printf("V bloku: %g\n", h_spectra_ref[512000].x);
/////////////////////////////////////////////////////////////////////////////
		float *h_real, *h_img; 
		float *d_coeff, *d_real, *d_img;
		Complex  *d_spectra;
		size_t free, total;
		CUresult res;

		h_real = new float[N];
		h_img = new float[N];

		for (int i = 0; i < N; i++){
		  h_real[i] = h_signal[i].x;
		  h_img[i] = h_signal[i].y;
		}
		if (debug) printf("x: %g y: %g\n", h_signal[0].x, h_signal[0].y);


		res = cuMemGetInfo(&free,&total);
		if (res != CUDA_SUCCESS) {printf("Memory read problem\n");}
		else if (debug) printf("free Memory: %lu %lu\n", free/1024/1024, total/1024/1024);
		//checkCudaErrors(cudaMalloc((void **) &d_signal,  sizeof(Complex)*filesize));




		checkCudaErrors(cudaMalloc((void **) &d_spectra, sizeof(Complex)*filesize));
		//checkCudaErrors(cudaMalloc((void **) &d_coeff,   sizeof(float)*nChannels*nTaps));
		checkCudaErrors(cudaMalloc((void **) &d_real,    sizeof(float)*filesize));
		checkCudaErrors(cudaMalloc((void **) &d_img,     sizeof(float)*filesize));
		checkCudaErrors(cudaMemcpyToSymbol(C_coeff,  h_coeff,   sizeof(float)*nChannels*nTaps));

		checkCudaErrors(cudaMemset(d_spectra,  0.0, sizeof(Complex)*filesize));

		//		checkCudaErrors(cudaMemcpy( d_coeff,   h_coeff,    nChannels*nTaps*sizeof(float), cudaMemcpyHostToDevice));
		//		checkCudaErrors(cudaMemcpy(d_signal,  h_signal,         filesize*sizeof(Complex), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_real,  h_real, filesize*sizeof(float), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_img,   h_img,  filesize*sizeof(float), cudaMemcpyHostToDevice));
		res = cuMemGetInfo(&free,&total);
		if (res != CUDA_SUCCESS) {printf("Memory read problem\n");}
		else if (debug) printf("free Memory: %lu %lu\n", free/1024/1024, total/1024/1024);


		//	    		complex2real<<<nBlocks,nChannels>>>(d_signal, d_real, d_img);
			//-----------------------------------------------------------------------
			//-----------------> launch the kernel
			dim3 gridSize(nBlocks - nTaps + 1, NUM_BLOCKS, 1);
			dim3 blockSize(nChannels/gridSize.y, 1, 1);
			int n_test = 1;
			float sum_time = 0.0f;
			if (TEST == 0) n_test = 1000;
			for (int i = 0; i < n_test; i++){
			  timer.Start();
			  Fir<<<gridSize, blockSize>>>(d_real, d_img, nTaps, nChannels, d_spectra);
			  timer.Stop();
			  sum_time += timer.Elapsed();
			}
			cudaDeviceSynchronize();
			checkCudaErrors(cudaGetLastError());
			printf("%d %g %i\n", nBlocks, sum_time /n_test, nChannels/gridSize.y);
			checkCudaErrors(cudaMemcpy(h_spectra,  d_spectra,  filesize*sizeof(Complex), cudaMemcpyDeviceToHost));
			//reference routine call - in utils_reference.h
			if (debug){
			  error = reference_code(h_spectra_ref, h_spectra, nChannels, nTaps);
				printf( "error = %g\n", error);
			}
			/*
		       	checkCudaErrors(cudaMemset(d_spectra,  0.0, sizeof(Complex)*filesize));
		       	timer.Start();
			Fir_01<<<gridSize, blockSize, 2*nTaps*nChannels*sizeof(float)>>>(d_real, d_img, d_coeff, nTaps, nChannels, d_spectra);
			timer.Stop();
			printf("%d %g %i\n", nBlocks, timer.Elapsed(), nChannels/gridSize.y);
			checkCudaErrors(cudaMemcpy(h_spectra,  d_spectra,  nChannels*sizeof(Complex), cudaMemcpyDeviceToHost));
			//reference routine call - in utils_reference.h
			if (debug){
			  error = reference_code(h_spectra_ref, h_spectra, nChannels, nBlocks);
				printf( "error = %g\n", error);
			}
			*/
	/*		test_01<<<M, N>>>(d_data, d_odata);
			checkCudaErrors(cudaMemcpy(h_odata, d_odata, sizeof(float)*N, cudaMemcpyDeviceToHost));
			printf("test_01: %g\n", h_odata[32]);
			*/

			// force the printf()s to flush
			
			
			//if (debug) cout << "Number of bl: " << bl << endl;
			//char str[200];
			//sprintf(str,"outtest_w_DFT.dat");
			//Save_data_spectra(str, h_spectra, filesize, nTaps, nChannels);


////////////////////////////////////////////////////////////////////////////
//  cuFFT
/////////////////////////////////////////////////////////////////////////////
/*
		//Create fft Plan
		cufftHandle plan;
		cufftSafeCall(cufftPlan1d(&plan, nChannels, CUFFT_C2C, nBlocks));

		//execute plan and copy back to host
		cufftSafeCall(cufftExecC2C(plan, (cufftComplex *)d_spectra, (cufftComplex *)d_spectra, CUFFT_FORWARD));
		checkCudaErrors(cudaMemcpy(h_spectra, d_spectra, nBlocks*nChannels*sizeof(Complex), cudaMemcpyDeviceToHost));
		
	       
		//write data to file
		//char str[200];
		//sprintf(str,"outtest_w_DFT.dat");
		//Save_data_spectra(str, h_spectra, filesize, nTaps, nChannels);

	    //Destroy the cuFFT plan
	    cufftSafeCall(cufftDestroy(plan));

/////////////////////////////////////////////////////////////////////////////
//  FFT reference
/////////////////////////////////////////////////////////////////////////////
		fftwf_plan p;
		fftwf_complex *in, *out;

		in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChannels);
		out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nChannels);
		p = fftwf_plan_dft_1d(nChannels, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		if (!p){
			exit(101);
		}

		for(int bl=nTaps-1;bl<nBlocks;bl++){
			// I need better memory arrangement
			for(int f=0; f<nChannels; f++){
				in[f][0] = h_spectra_ref[bl*nChannels+f].x;
				in[f][1] = h_spectra_ref[bl*nChannels+f].y;
			}
			fftwf_execute(p);
			for(int f=0;f<nChannels;f++){
				h_spectra_ref[bl*nChannels+f].x = out[f][0];
				h_spectra_ref[bl*nChannels+f].y = out[f][1];
			}
		}
		
		//printf("%g\n",h_spectra_ref[257].x - h_spectra[0].x);
		//char str[200];
		//sprintf(str,"outtestCPU_w_DFT.dat");
		//Save_data_spectra(str, h_spectra_ref, filesize, nTaps, nChannels);
		//reference routine call - in utils_reference.h
		if (debug) {
		  error = reference_code(h_spectra_ref, h_spectra, nChannels, nBlocks);
			printf( "error = %g\n", error);
		}

	    fftwf_destroy_plan(p);
	    fftwf_free(in);
	    fftwf_free(out);
*/
/////////////////////////////////////////////////////////////////////////////
//  Clean-up the rest -- be a good kid!
/////////////////////////////////////////////////////////////////////////////

		delete[] h_coeff;
		delete[] h_spectra_ref;
		delete[] w_buffer;
		delete[] h_signal;
		delete[] h_spectra;
		delete[] h_real;
		delete[] h_img;
		//   checkCudaErrors(cudaFree(C_coeff));
	    //	    checkCudaErrors(cudaFree(d_signal));
       	    checkCudaErrors(cudaFree(d_spectra));
	    checkCudaErrors(cudaFree(d_real));
	    checkCudaErrors(cudaFree(d_img));

	    //	    cudaDeviceReset();
	    // say bye bye
	    if (debug) printf("That's all Folks!\n");
	    
	    //	    cudaDeviceReset();
	return 0;
}
