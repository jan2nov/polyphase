#ifndef FIR_WINDOW
#define FIR_WINDOW
#endif

extern int THREADS;
//TODO:
// From intrinsic implementation Mk1 I do not need to parallelize the code by blocks of a size nBlocks/nThreads. This was due to rotary buffer and data reuse.

//Fourth intrinsic implementation. I introduced unified data approach. This requaries new coeff coefficients whose are twice the size of normal coeff. The main advantage of this is I do not need to transfor data to and from real,imaginary format.
//Tested? YES
void FIR_para_INT_Mk3(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nBlocks, int nThreads);
void FIR_serial_INT_Mk3(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nBlocks);

void FIR_para_INT_Mk4(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nBlocks, int nThreads);
void FIR_serial_INT_Mk4(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nBlocks);
void FIR_serial_Mk1_INT_Mk4(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nBlocks);

void FIR_para_INT_Mk5(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nBlocks, int nThreads);
void FIR_para_INT_Mk6(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nBlocks,int nThreads);

void FIR_para_INT_Mk4_256(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nBlocks, int nThreads);
void FIR_serial_INT_Mk4_256(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nBlocks);


































