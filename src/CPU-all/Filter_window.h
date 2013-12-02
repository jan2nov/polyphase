#ifndef FIR_WINDOW
#define FIR_WINDOW
#endif

/*
#include <math.h>
#include <fftw3.h>
const double M_PI=3.141592654;
*/

double besselI0(double x);


void W_kaiser(int n, double beta, double *window);


// Guassian window function
void W_gaussian(int n, double a, double *window);


void W_hamming(unsigned n, double *window);


void W_blackman(unsigned n, double* d);


bool interpolate(const double* x, const double* y, unsigned int nX, unsigned int n, double* result);


unsigned int nextPowerOf2(unsigned int n);


void Generate_FIR_Filter(unsigned int n, double w, const double* window, double* result);


bool Window_coefficients(unsigned int nTaps, unsigned int nChannels, int windowtype, float *coeff);

