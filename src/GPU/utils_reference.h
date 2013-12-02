float reference_code(float2 *h_spectra_ref, float2 *h_spectra, int nChannels, unsigned int nBlocks){

	float diff, error_norm = 0;
	int nTaps = 8;

	for (int j = 0; j < nBlocks - nTaps + 1; j++){
	  for (int i = 0; i < nChannels; i++) {
	    diff = h_spectra_ref[i + 7*nChannels + j*nChannels].x - h_spectra[j*nChannels + i].x;
	      error_norm += diff * diff;
	      //printf("Ahoj %g %g\n", h_spectra_ref[i + 7*nChannels].x, h_spectra[i].x);
	      diff = h_spectra_ref[i + 7*nChannels + j*nChannels].y - h_spectra[j*nChannels + i].y;
	      error_norm += diff * diff;
	  }
	}
	error_norm = (float)sqrt((double)error_norm);

	return error_norm;
}
