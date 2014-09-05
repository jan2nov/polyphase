__global__ void Fir_SpB_shared(float2* __restrict__  d_data, float* __restrict__ d_coeff, int nTaps, int nChannels, int yshift, float2* __restrict__ d_spectra) {
	int t = 0;
	int bl= 4*blockIdx.x*nChannels;
	int ypos = blockDim.x*blockIdx.y + yshift;
	float2 ftemp1 = make_float2(0.0,0.0);
	float2 ftemp2 = make_float2(0.0,0.0);
	float2 ftemp3 = make_float2(0.0,0.0);
	float2 ftemp4 = make_float2(0.0,0.0);
	float temp;

	for(t=ypos + threadIdx.x;t<(nTaps)*nChannels;t+=nChannels){
		temp = d_coeff[t]; 
		ftemp1.x = __fmaf_rn(temp,d_data[bl+t].x,ftemp1.x);
		ftemp1.y = __fmaf_rn(temp,d_data[bl+t].y,ftemp1.y);
		ftemp2.x = __fmaf_rn(temp,d_data[bl+nChannels+t].x,ftemp2.x);
		ftemp2.y = __fmaf_rn(temp,d_data[bl+nChannels+t].y,ftemp2.y);
		ftemp3.x = __fmaf_rn(temp,d_data[bl+2*nChannels+t].x,ftemp3.x);
		ftemp3.y = __fmaf_rn(temp,d_data[bl+2*nChannels+t].y,ftemp3.y);
		ftemp4.x = __fmaf_rn(temp,d_data[bl+3*nChannels+t].x,ftemp4.x);
		ftemp4.y = __fmaf_rn(temp,d_data[bl+3*nChannels+t].y,ftemp4.y);
	}

	t=bl + ypos + threadIdx.x;
	d_spectra[t]=ftemp1;
	d_spectra[t+nChannels]=ftemp2;
	d_spectra[t+2*nChannels]=ftemp3;
	d_spectra[t+3*nChannels]=ftemp4;
	return;
}
