#include <cstdlib>
#include <math.h>
#include <fftw3.h>

const double M_PI=3.141592654;


double besselI0(double x) {
    // Parameters of the polynomial approximation.
    const double p1 = 1.0, p2 = 3.5156229, p3 = 3.0899424, p4 = 1.2067492,
            p5 = 0.2659732, p6 = 3.60768e-2,  p7 = 4.5813e-3;

    const double q1 = 0.39894228, q2 = 1.328592e-2, q3 = 2.25319e-3,
            q4 = -1.57565e-3, q5 = 9.16281e-3, q6 = -2.057706e-2,
            q7 = 2.635537e-2, q8 = -1.647633e-2, q9 = 3.92377e-3;

    const double k1 = 3.75;

    double ax = abs(x);
    double y = 0, result = 0;

    if (ax < k1) {
        double xx = x / k1;
        y = xx * xx;
        result = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))));
    }
    else {
        y = k1 / ax;
        result = (exp(ax)/sqrt(ax))*
                (q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
    }
    return result;
}


void W_kaiser(int n, double beta, double *window)
{
    if (n == 1) {
        window[0] = 1.0;
        return;
    }

    int m = n - 1;

    for (int i = 0; i < n; i++) {
        double k = 2.0 * beta / m * sqrt(double(i * (m - i)));
        window[i] = besselI0(k) / besselI0(beta);
    }
}


// Guassian window function
void W_gaussian(int n, double a, double *window)
{
    int index = 0;
    for (int i=-(n-1); i<=n-1; i+=2) {
        window[index++] = exp( -0.5 * pow(( a/n * i), 2) );
    }
}


void W_hamming(unsigned n, double *window)
{
    if (n == 1) {
        window[0] = 1.0;
        return;
    }
    unsigned m = n - 1;
    for(unsigned i = 0; i < n; i++) {
        window[i] = 0.54 - 0.46 * cos((2.0 * M_PI * i) / m);
    }
}




void W_blackman(unsigned n, double* d)
{
    if(n == 1) {
        d[0] = 1.0;
        return;
    }
    unsigned m = n - 1;
    for(unsigned i = 0; i < n; i++) {
        double k = i / m;
        d[i] = 0.42 - 0.5 * cos(2.0 * M_PI * k) + 0.08 * cos(4.0 * M_PI * k);
    }
}


bool interpolate(const double* x, const double* y, unsigned int nX, unsigned int n, double* result) {
    unsigned int nextX = 0;
    unsigned int index = 0;

    for(double interpolatedX = 0.0; interpolatedX <= 1.0; interpolatedX += 1.0/(n-1), index++) {

        while(x[nextX] <= interpolatedX && nextX < nX - 1) {
            nextX++;
        }

        if(nextX == 0) {
            return(1);
        }

        double prevXVal = x[nextX-1];
        double nextXVal = x[nextX];
        double prevYVal = y[nextX-1];
        double nextYVal = y[nextX];

        double rc = (nextYVal - prevYVal) / (nextXVal - prevXVal);
        double newVal = prevYVal + (interpolatedX - prevXVal) * rc;

        result[index] = newVal;
    }
	return(0);
}



unsigned int nextPowerOf2(unsigned int n) {
    unsigned int res = 1;
    while(true) {
        if(res >= n) { return res; }
        res *= 2;
    }
}

void Generate_FIR_Filter(unsigned int n, double w, const double* window, double* result) {
	bool error=0;
    // make sure grid is big enough for the window
    // the grid must be at least (n+1)/2
    // for all filters where the order is a power of two minus 1, grid_n = n+1;
    unsigned int grid_n = nextPowerOf2(n + 1);
    unsigned int ramp_n = 2; // grid_n / 20;

    // Apply ramps to discontinuities
    // this is a low pass filter
    // maybe we can omit the "w, 0" point?
    // I did observe a small difference
    double f[] = {0.0, w-ramp_n/grid_n/2.0, w, w+ramp_n/grid_n/2.0, 1.0};
    double m[] = {1.0, 1.0, 0.0, 0.0, 0.0};

    // grid is a 1-D array with grid_n+1 points. Values are 1 in filter passband, 0 otherwise
    double* grid = (double*) malloc((grid_n + 1) * sizeof(double));

    // interpolate between grid points
    error=interpolate(f, m, 5 /* length of f and m arrays */ , grid_n+1, grid);
	if (error) {
		printf("Error in interpolation!");
		exit(10);
	}

    // the grid we do an ifft on is:
    // grid appended with grid_n*2 zeros
    // appended with original grid values from indices grid_n..2, i.e., the values in reverse order
    // (note, arrays start at 1 in octave!)
    // the input for the ifft is of size 4*grid_n
    // input = [grid ; zeros(grid_n*2,1) ;grid(grid_n:-1:2)];

    fftwf_complex* cinput  = (fftwf_complex*) fftwf_malloc(grid_n*4*sizeof(fftwf_complex));
    fftwf_complex* coutput = (fftwf_complex*) fftwf_malloc(grid_n*4*sizeof(fftwf_complex));

    if(cinput == NULL || coutput == NULL) {
        printf("Cannot allocate buffers FFTW");
		exit(11);
    }

    // wipe imaginary part
    for(unsigned i=0; i<grid_n*4; i++) {
        cinput[i][1] = 0.0;
    }

    // copy first part of grid
    for(unsigned i=0; i<grid_n+1; i++) {
        cinput[i][0] = float(grid[i]);
    }

    // append zeros
    for(unsigned i=grid_n+1; i<=grid_n*3; i++) {
        cinput[i][0] = 0.0;
    }

    // now append the grid in reverse order
    for(unsigned i=grid_n-1, index=0; i >=1; i--, index++) {
        cinput[grid_n*3+1 + index][0] = float(grid[i]);
    }

    fftwf_plan plan = fftwf_plan_dft_1d(grid_n*4, cinput, coutput,
            FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);

    unsigned index = 0;
    for(unsigned i=4*grid_n-n; i<4*grid_n; i+=2) {
        result[index] = coutput[i][0];
        index++;
    }

    for(unsigned i=1; i<=n; i+=2) {
        result[index] = coutput[i][0];
        index++;
    }

    fftwf_destroy_plan(plan);
    fftwf_free(cinput);
    fftwf_free(coutput);

    // multiply with window
    for(unsigned i=0; i<=n; i++) {
        result[i] *= window[i];
    }

    // normalize
    double factor = result[n/2];
    for(unsigned i=0; i<=n; i++) {
        result[i] /= factor;
    }

    free(grid);
}



bool Window_coefficients(unsigned int nTaps, unsigned int nChannels, int windowtype, float *coeff) {
    unsigned n = nChannels * nTaps;
    double* window = new double[n];

	// Error checking
	if(windowtype < 1 || windowtype > 4) {
		printf("Wrong windowtype");
		return(1);
	}

	if (windowtype == 1) { //kaiser
		double beta = 9.0695;
		W_kaiser(n, beta, window);
	}
	else if (windowtype == 2) {//gaussian
		double alpha = 3.5;
		W_gaussian(n, alpha, window);
	}
	else if (windowtype == 3) {//blackman
		W_blackman(n, window);
	}
	else if (windowtype == 4) {//hamming
		W_hamming(n, window);
	}

    double* result = new double[n];
    Generate_FIR_Filter(n - 1, 1.0 / nChannels, window, result);

    delete[] window;

    for(unsigned int t = 0; t < nTaps; ++t) {
        for(unsigned int c = 0; c < nChannels; ++c) {
            unsigned int index = t * nChannels + c;
            coeff[index] = (float) (result[index] / nChannels);
            if (c%2 == 0)
                coeff[index] = (float) (result[index] / nChannels);
            else
                coeff[index] = (float) (-result[index] / nChannels);
        }
    }

    delete[] result;
	return(0);
}


