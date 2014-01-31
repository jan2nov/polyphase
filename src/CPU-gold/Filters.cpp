#include <stdio.h>
#include <omp.h>
#include <xmmintrin.h>
#include <immintrin.h>
#include <string.h>
// open mp is defined by OMP_H
// _XMMINTRIN_H_INCLUDED
int THREADS;



void Setup_buffer(float *w_buffer, float *real, float *img, unsigned int *oldesttap, int nChannels, int nTaps){
	*oldesttap=7;

	for(int i=0;i<nTaps;i++){
		memcpy(&w_buffer[i*nChannels*2],&real[i*nChannels],sizeof(float)*nChannels);
		memcpy(&w_buffer[i*nChannels*2+nChannels],&img[i*nChannels],sizeof(float)*nChannels);
	}

}

void Update_buffer(float *w_buffer, float *real, float *img, unsigned int *oldesttap, int nChannels, int nTaps){
	unsigned int itemp=*oldesttap;
	memcpy(&w_buffer[itemp*nChannels*2],real,sizeof(float)*nChannels);
	memcpy(&w_buffer[itemp*nChannels*2+nChannels],img,sizeof(float)*nChannels);
	
	itemp=(itemp+1)%nTaps;
	*oldesttap=itemp;
}


//**************
// SERIAL CODE
//**************
void FIR_serial(float *real, float *img, float *spectra_real, float *spectra_img, float *coeff, unsigned int nTaps, unsigned int nChannels, unsigned int nBlocks){
	float *w_buffer;
	unsigned int tap=0;
	unsigned int oldesttap;
	
	w_buffer=new float[nChannels*nTaps*2];

	Setup_buffer(w_buffer,real,img,&oldesttap,nChannels,nTaps);

	for(int bl=nTaps-1;bl<nBlocks;bl++){
		
		Update_buffer(w_buffer,&real[bl*nChannels],&img[bl*nChannels],&oldesttap,nChannels,nTaps);

		for(int t=0;t<nTaps;t++){
			tap=(oldesttap+t)%nTaps;
			for(int c=0;c<nChannels;c++){
				spectra_real[bl*nChannels+c]+=coeff[t*nChannels+c]*w_buffer[tap*2*nChannels+c];
				spectra_img[bl*nChannels+c]+=coeff[t*nChannels+c]*w_buffer[tap*2*nChannels+c+nChannels];
			}//nChannels
		}//nTaps
	} //nBlocks
}


//**************
// SERIAL CODE WITH INTRINSIC Mk3 for architecture with register size = 128bit
//**************
void FIR_serial_INT_Mk3(float *input_data, float *spectra, float *coeff, unsigned int nTaps, unsigned int nChannels, unsigned int nBlocks){
	//number of channels must by dividable by factor of 8
	int c,t,bl,tap,th_id,nthreads,block_step;
	int itemp1,itemp2;
	nChannels=nChannels*2;
	
	//registers
	int FpR=4,RpCl=4;
	int REG=nChannels/(RpCl*FpR);
	__m128 i_data[4];// because of the size of the cacheline 64byte
	__m128 i_coeff[4];
	__m128 i_spectra[4];
	__m128 i_temp[4];
	
		for(bl=0;bl<nBlocks-(nTaps-1);bl++){
			// cycle through channels in multiples of 16
			for(c=0;c<REG;c++){
				// Zeroing!
				i_spectra[0]=_mm_setzero_ps();i_spectra[1]=_mm_setzero_ps();i_spectra[2]=_mm_setzero_ps();i_spectra[3]=_mm_setzero_ps();
				
				for(t=0;t<nTaps;t++){
					//Loading data
					//_mm_loadu_ps is used, this is before we implemented alignment of memory.  Intrinsic _mm_loadu_ps do not assume alignment of memory
					i_coeff[0]=_mm_loadu_ps(&coeff[(RpCl*c)*FpR+t*nChannels]);
					i_coeff[1]=_mm_loadu_ps(&coeff[(RpCl*c+1)*FpR+t*nChannels]);
					i_coeff[2]=_mm_loadu_ps(&coeff[(RpCl*c+2)*FpR+t*nChannels]);
					i_coeff[3]=_mm_loadu_ps(&coeff[(RpCl*c+3)*FpR+t*nChannels]);
										
					i_data[0]=_mm_loadu_ps(&input_data[(t+bl)*nChannels+(RpCl*c)*FpR]);
					i_data[1]=_mm_loadu_ps(&input_data[(t+bl)*nChannels+(RpCl*c+1)*FpR]);
					i_data[2]=_mm_loadu_ps(&input_data[(t+bl)*nChannels+(RpCl*c+2)*FpR]);
					i_data[3]=_mm_loadu_ps(&input_data[(t+bl)*nChannels+(RpCl*c+3)*FpR]);
					
					//Computation
					i_temp[0]=_mm_mul_ps(i_data[0],i_coeff[0]);
					i_temp[1]=_mm_mul_ps(i_data[1],i_coeff[1]);
					i_temp[2]=_mm_mul_ps(i_data[2],i_coeff[2]);
					i_temp[3]=_mm_mul_ps(i_data[3],i_coeff[3]);
					
					i_spectra[0]=_mm_add_ps(i_temp[0],i_spectra[0]);
					i_spectra[1]=_mm_add_ps(i_temp[1],i_spectra[1]);
					i_spectra[2]=_mm_add_ps(i_temp[2],i_spectra[2]);
					i_spectra[3]=_mm_add_ps(i_temp[3],i_spectra[3]);			
				} // for nTaps
				// here we are making use of the fact that we do not need the resulting spectra and we are bypassing cache and writing directly into memory. This should avoid unnecessary cache writing and related data un-caching.
				_mm_stream_ps(&spectra[(RpCl*c)*FpR + bl*nChannels],i_spectra[0]);
				_mm_stream_ps(&spectra[(RpCl*c+1)*FpR + bl*nChannels],i_spectra[1]);
				_mm_stream_ps(&spectra[(RpCl*c+2)*FpR + bl*nChannels],i_spectra[2]);
				_mm_stream_ps(&spectra[(RpCl*c+3)*FpR + bl*nChannels],i_spectra[3]);
			}// nChannels
		}// for nBlocks
}


//**************
// PARALLEL CODE WITH INTRINSIC Mk3 for architecture with register size = 128bit
//**************
void FIR_para_INT_Mk3(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nBlocks, int nThreads){
	//number of channels must by dividable by factor of 8
	int c,t,bl,th_id,block_step;
	nChannels=nChannels*2;
	
	//registers
	int FpR=4,RpCl=4;
	int REG=nChannels/(RpCl*FpR);
	__m128 i_data[4];// because of the size of the cacheline 64byte
	__m128 i_coeff[4];
	__m128 i_spectra[4];
	__m128 i_temp[4];
	
	omp_set_num_threads(nThreads);
	#pragma omp parallel shared(input_data,spectra,coeff) private(bl,t,c,i_data,i_coeff,i_spectra,i_temp)
	{
		//th_id = omp_get_thread_num();
		//nThreads = omp_get_num_threads();
		block_step=nBlocks/nThreads;
		
		#pragma omp for schedule(static,block_step)
		for(bl=0;bl<nBlocks-(nTaps-1);bl++){
			// cycle through channels in multiples of 16
			for(c=0;c<REG;c++){
				// Zeroing!
				i_spectra[0]=_mm_setzero_ps();i_spectra[1]=_mm_setzero_ps();i_spectra[2]=_mm_setzero_ps();i_spectra[3]=_mm_setzero_ps();
				
				for(t=0;t<nTaps;t++){
					//Loading data
					//_mm_loadu_ps is used, this is before we implemented alignment of memory.  Intrinsic _mm_loadu_ps do not assume alignment of memory					
					i_coeff[0]=_mm_loadu_ps(&coeff[(RpCl*c)*FpR+t*nChannels]);
					i_coeff[1]=_mm_loadu_ps(&coeff[(RpCl*c+1)*FpR+t*nChannels]);
					i_coeff[2]=_mm_loadu_ps(&coeff[(RpCl*c+2)*FpR+t*nChannels]);
					i_coeff[3]=_mm_loadu_ps(&coeff[(RpCl*c+3)*FpR+t*nChannels]);
										
					i_data[0]=_mm_loadu_ps(&input_data[(t+bl)*nChannels+(RpCl*c)*FpR]);
					i_data[1]=_mm_loadu_ps(&input_data[(t+bl)*nChannels+(RpCl*c+1)*FpR]);
					i_data[2]=_mm_loadu_ps(&input_data[(t+bl)*nChannels+(RpCl*c+2)*FpR]);
					i_data[3]=_mm_loadu_ps(&input_data[(t+bl)*nChannels+(RpCl*c+3)*FpR]);
					
					i_temp[0]=_mm_mul_ps(i_data[0],i_coeff[0]);
					i_temp[1]=_mm_mul_ps(i_data[1],i_coeff[1]);
					i_temp[2]=_mm_mul_ps(i_data[2],i_coeff[2]);
					i_temp[3]=_mm_mul_ps(i_data[3],i_coeff[3]);
					
					i_spectra[0]=_mm_add_ps(i_temp[0],i_spectra[0]);
					i_spectra[1]=_mm_add_ps(i_temp[1],i_spectra[1]);
					i_spectra[2]=_mm_add_ps(i_temp[2],i_spectra[2]);
					i_spectra[3]=_mm_add_ps(i_temp[3],i_spectra[3]);			
				} // for nTaps
				// here we are making use of the fact that we do not need the resulting spectra and we are bypassing cache and writing directly into memory. This should avoid unnecessary cache writing and related data un-caching.
				_mm_stream_ps(&spectra[(RpCl*c)*FpR + bl*nChannels],i_spectra[0]);
				_mm_stream_ps(&spectra[(RpCl*c+1)*FpR + bl*nChannels],i_spectra[1]);
				_mm_stream_ps(&spectra[(RpCl*c+2)*FpR + bl*nChannels],i_spectra[2]);
				_mm_stream_ps(&spectra[(RpCl*c+3)*FpR + bl*nChannels],i_spectra[3]);
			}// nChannels
		}// for nBlocks
		
	}// parallel block	
}

//**************
// PARALLEL CODE WITH INTRINSIC Mk4 for architecture with register size = 128bit
//**************
void FIR_para_INT_Mk4(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nBlocks, int nThreads){
	//number of channels must by dividable by factor of 8
	int c,t,bl,th_id,block_step;
	nChannels=nChannels*2;
	
	//registers
	int FpR=4,RpCl=4;
	int REG=nChannels/(RpCl*FpR);
	__m128 i_data[4];// because of the size of the cacheline 64byte
	__m128 i_coeff[4];
	__m128 i_spectra[4];
	__m128 i_temp[4];
	
	omp_set_num_threads(nThreads);
	#pragma omp parallel shared(input_data,spectra,coeff) private(bl,t,c,i_data,i_coeff,i_spectra,i_temp)
	{
		//th_id = omp_get_thread_num();
		//nThreads = omp_get_num_threads();
		block_step=nBlocks/nThreads;
		
		#pragma omp for schedule(static,block_step)
		for(bl=0;bl<nBlocks-(nTaps-1);bl++){
			for(c=0;c<REG;c++){
				// Zeroing!
				i_spectra[0]=_mm_setzero_ps();i_spectra[1]=_mm_setzero_ps();i_spectra[2]=_mm_setzero_ps();i_spectra[3]=_mm_setzero_ps();
				
				for(t=0;t<nTaps;t++){
					//Loading data
					// This is the only change in algorithm with regards to Mk3 function. We are using  _mm_load_ps which requires aligned memory.
					i_coeff[0]=_mm_load_ps(&coeff[(RpCl*c)*FpR+t*nChannels]);
					i_coeff[1]=_mm_load_ps(&coeff[(RpCl*c+1)*FpR+t*nChannels]);
					i_coeff[2]=_mm_load_ps(&coeff[(RpCl*c+2)*FpR+t*nChannels]);
					i_coeff[3]=_mm_load_ps(&coeff[(RpCl*c+3)*FpR+t*nChannels]);
										
					i_data[0]=_mm_load_ps(&input_data[(t+bl)*nChannels+(RpCl*c)*FpR]);
					i_data[1]=_mm_load_ps(&input_data[(t+bl)*nChannels+(RpCl*c+1)*FpR]);
					i_data[2]=_mm_load_ps(&input_data[(t+bl)*nChannels+(RpCl*c+2)*FpR]);
					i_data[3]=_mm_load_ps(&input_data[(t+bl)*nChannels+(RpCl*c+3)*FpR]);
					
					i_temp[0]=_mm_mul_ps(i_data[0],i_coeff[0]);
					i_temp[1]=_mm_mul_ps(i_data[1],i_coeff[1]);
					i_temp[2]=_mm_mul_ps(i_data[2],i_coeff[2]);
					i_temp[3]=_mm_mul_ps(i_data[3],i_coeff[3]);
					
					i_spectra[0]=_mm_add_ps(i_temp[0],i_spectra[0]);
					i_spectra[1]=_mm_add_ps(i_temp[1],i_spectra[1]);
					i_spectra[2]=_mm_add_ps(i_temp[2],i_spectra[2]);
					i_spectra[3]=_mm_add_ps(i_temp[3],i_spectra[3]);			
				} // for nTaps
				_mm_stream_ps(&spectra[(RpCl*c)*FpR + bl*nChannels],i_spectra[0]);
				_mm_stream_ps(&spectra[(RpCl*c+1)*FpR + bl*nChannels],i_spectra[1]);
				_mm_stream_ps(&spectra[(RpCl*c+2)*FpR + bl*nChannels],i_spectra[2]);
				_mm_stream_ps(&spectra[(RpCl*c+3)*FpR + bl*nChannels],i_spectra[3]);
			}// nChannels
		}// for nBlocks
		
	}// parallel block	
}

//**************
// SERIAL CODE WITH INTRINSIC Mk4 for architecture with register size = 128bit
//**************
void FIR_serial_INT_Mk4(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nBlocks){
	//number of channels must by dividable by factor of 8
	int c,t,bl,th_id,block_step;
	nChannels=nChannels*2;
	
	//registers
	int FpR=4,RpCl=4;
	int REG=nChannels/(RpCl*FpR);
	__m128 i_data[4];// because of the size of the cacheline 64byte
	__m128 i_coeff[4];
	__m128 i_spectra[4];
	__m128 i_temp[4];
	
		for(bl=0;bl<nBlocks-(nTaps-1);bl++){
			for(c=0;c<REG;c++){
				// Zeroing!
				i_spectra[0]=_mm_setzero_ps();i_spectra[1]=_mm_setzero_ps();i_spectra[2]=_mm_setzero_ps();i_spectra[3]=_mm_setzero_ps();
				
				for(t=0;t<nTaps;t++){
					//Loading data
					
					i_coeff[0]=_mm_load_ps(&coeff[(RpCl*c)*FpR+t*nChannels]);
					i_coeff[1]=_mm_load_ps(&coeff[(RpCl*c+1)*FpR+t*nChannels]);
					i_coeff[2]=_mm_load_ps(&coeff[(RpCl*c+2)*FpR+t*nChannels]);
					i_coeff[3]=_mm_load_ps(&coeff[(RpCl*c+3)*FpR+t*nChannels]);
										
					i_data[0]=_mm_load_ps(&input_data[(t+bl)*nChannels+(RpCl*c)*FpR]);
					i_data[1]=_mm_load_ps(&input_data[(t+bl)*nChannels+(RpCl*c+1)*FpR]);
					i_data[2]=_mm_load_ps(&input_data[(t+bl)*nChannels+(RpCl*c+2)*FpR]);
					i_data[3]=_mm_load_ps(&input_data[(t+bl)*nChannels+(RpCl*c+3)*FpR]);
					
					i_temp[0]=_mm_mul_ps(i_data[0],i_coeff[0]);
					i_temp[1]=_mm_mul_ps(i_data[1],i_coeff[1]);
					i_temp[2]=_mm_mul_ps(i_data[2],i_coeff[2]);
					i_temp[3]=_mm_mul_ps(i_data[3],i_coeff[3]);
					
					i_spectra[0]=_mm_add_ps(i_temp[0],i_spectra[0]);
					i_spectra[1]=_mm_add_ps(i_temp[1],i_spectra[1]);
					i_spectra[2]=_mm_add_ps(i_temp[2],i_spectra[2]);
					i_spectra[3]=_mm_add_ps(i_temp[3],i_spectra[3]);			
				} // for nTaps
				_mm_stream_ps(&spectra[(RpCl*c)*FpR + bl*nChannels],i_spectra[0]);
				_mm_stream_ps(&spectra[(RpCl*c+1)*FpR + bl*nChannels],i_spectra[1]);
				_mm_stream_ps(&spectra[(RpCl*c+2)*FpR + bl*nChannels],i_spectra[2]);
				_mm_stream_ps(&spectra[(RpCl*c+3)*FpR + bl*nChannels],i_spectra[3]);
			}// nChannels
		}// for nBlocks
}


//**************
// SERIAL CODE WITH INTRINSIC Mk4 for architecture with register size = 128bit
//**************
// This function is a attempt to optimize the code fore cache utilization of size 32kB. This is why the calculation is divided into chunks of 32 iteration assuming 512 channels and 8 taps
void FIR_serial_Mk1_INT_Mk4(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nBlocks){
	//number of channels must by dividable by factor of 256
	int i,c,t,bl,block_step;
	int itemp;
	nChannels=nChannels*2;
	
	//registers
	int FpR=4,RpCl=4;
	int REG=nChannels/(RpCl*FpR);
	__m128 i_data[4];// because of the size of the cacheline 64byte
	__m128 i_coeff[4];
	__m128 i_spectra[4];
	__m128 i_temp[4];
	
	int INNER=32;
	int OUTER=REG/INNER;
	for(i=0;i<OUTER;i++){
		for(bl=0;bl<nBlocks-(nTaps-1);bl++){
			for(c=0;c<INNER;c++){
				// Zeroing!
				i_spectra[0]=_mm_setzero_ps();i_spectra[1]=_mm_setzero_ps();i_spectra[2]=_mm_setzero_ps();i_spectra[3]=_mm_setzero_ps();
				
				for(t=0;t<nTaps;t++){
					//Loading data
					i_coeff[0]=_mm_load_ps(&coeff[(RpCl*c)*FpR+t*nChannels + i*INNER*RpCl*FpR]);
					i_coeff[1]=_mm_load_ps(&coeff[(RpCl*c+1)*FpR+t*nChannels + i*INNER*RpCl*FpR]);
					i_coeff[2]=_mm_load_ps(&coeff[(RpCl*c+2)*FpR+t*nChannels + i*INNER*RpCl*FpR]);
					i_coeff[3]=_mm_load_ps(&coeff[(RpCl*c+3)*FpR+t*nChannels + i*INNER*RpCl*FpR]);
										
					i_data[0]=_mm_load_ps(&input_data[(bl+t)*nChannels + (RpCl*c)*FpR + i*INNER*RpCl*FpR]);
					i_data[1]=_mm_load_ps(&input_data[(bl+t)*nChannels + (RpCl*c+1)*FpR + i*INNER*RpCl*FpR]);
					i_data[2]=_mm_load_ps(&input_data[(bl+t)*nChannels + (RpCl*c+2)*FpR + i*INNER*RpCl*FpR]);
					i_data[3]=_mm_load_ps(&input_data[(bl+t)*nChannels + (RpCl*c+3)*FpR + i*INNER*RpCl*FpR]);
					
					i_temp[0]=_mm_mul_ps(i_data[0],i_coeff[0]);
					i_temp[1]=_mm_mul_ps(i_data[1],i_coeff[1]);
					i_temp[2]=_mm_mul_ps(i_data[2],i_coeff[2]);
					i_temp[3]=_mm_mul_ps(i_data[3],i_coeff[3]);
					
					i_spectra[0]=_mm_add_ps(i_temp[0],i_spectra[0]);
					i_spectra[1]=_mm_add_ps(i_temp[1],i_spectra[1]);
					i_spectra[2]=_mm_add_ps(i_temp[2],i_spectra[2]);
					i_spectra[3]=_mm_add_ps(i_temp[3],i_spectra[3]);			
				} // for nTaps
				_mm_stream_ps(&spectra[(RpCl*c)*FpR + i*INNER*RpCl*FpR + bl*nChannels],i_spectra[0]);
				_mm_stream_ps(&spectra[(RpCl*c+1)*FpR + i*INNER*RpCl*FpR + bl*nChannels],i_spectra[1]);
				_mm_stream_ps(&spectra[(RpCl*c+2)*FpR + i*INNER*RpCl*FpR + bl*nChannels],i_spectra[2]);
				_mm_stream_ps(&spectra[(RpCl*c+3)*FpR + i*INNER*RpCl*FpR + bl*nChannels],i_spectra[3]);
			}//inner
		}// for nBlocks
	}// nChannels
}

//**************
// PARALLEL CODE WITH INTRINSIC Mk4 BASED ON SERIAL Mk1 for architecture with register size = 128bit
//**************
//parallel version of the cache optimized code
void FIR_para_INT_Mk6(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nBlocks,int nThreads){
	//number of channels must by dividable by factor of 256
	int i,c,t,bl,th_id,block_step;
	int start,end;
	nChannels=nChannels*2;
	
	//registers
	int FpR=4,RpCl=4;
	int REG=nChannels/(RpCl*FpR);
	__m128 i_data[4];// because of the size of the cacheline 64byte
	__m128 i_coeff[4];
	__m128 i_spectra[4];
	__m128 i_temp[4];
	
	// Outer for
	// Inner for
	int INNER=32;
	int OUTER=REG/INNER;
	omp_set_num_threads(nThreads);
	#pragma omp parallel shared(input_data,spectra,coeff) private(start,end,i,bl,t,c,i_data,i_coeff,i_spectra,i_temp)
	{
		th_id = omp_get_thread_num();
		block_step=nBlocks/nThreads;
		if (block_step==0){
			if ( th_id<(nBlocks-(nTaps-1)) ) {
				start=th_id;
				end=th_id+1;
			}
		}
		else {
			end=(th_id+1)*block_step;
			if (th_id==nThreads-1) end=end-(nTaps-1);
			start=th_id*block_step;
		}
		
		for(i=0;i<OUTER;i++){
		
			for(bl=start;bl<end;bl++){
				for(c=0;c<INNER;c++){
					// Zeroing!
					i_spectra[0]=_mm_setzero_ps();i_spectra[1]=_mm_setzero_ps();i_spectra[2]=_mm_setzero_ps();i_spectra[3]=_mm_setzero_ps();
					
					for(t=0;t<nTaps;t++){
						//Loading data
						i_coeff[0]=_mm_load_ps(&coeff[(RpCl*c)*FpR+t*nChannels + i*INNER*RpCl*FpR]);
						i_coeff[1]=_mm_load_ps(&coeff[(RpCl*c+1)*FpR+t*nChannels + i*INNER*RpCl*FpR]);
						i_coeff[2]=_mm_load_ps(&coeff[(RpCl*c+2)*FpR+t*nChannels + i*INNER*RpCl*FpR]);
						i_coeff[3]=_mm_load_ps(&coeff[(RpCl*c+3)*FpR+t*nChannels + i*INNER*RpCl*FpR]);
											
						i_data[0]=_mm_load_ps(&input_data[(bl+t)*nChannels + (RpCl*c)*FpR + i*INNER*RpCl*FpR]);
						i_data[1]=_mm_load_ps(&input_data[(bl+t)*nChannels + (RpCl*c+1)*FpR + i*INNER*RpCl*FpR]);
						i_data[2]=_mm_load_ps(&input_data[(bl+t)*nChannels + (RpCl*c+2)*FpR + i*INNER*RpCl*FpR]);
						i_data[3]=_mm_load_ps(&input_data[(bl+t)*nChannels + (RpCl*c+3)*FpR + i*INNER*RpCl*FpR]);
						
						i_temp[0]=_mm_mul_ps(i_data[0],i_coeff[0]);
						i_temp[1]=_mm_mul_ps(i_data[1],i_coeff[1]);
						i_temp[2]=_mm_mul_ps(i_data[2],i_coeff[2]);
						i_temp[3]=_mm_mul_ps(i_data[3],i_coeff[3]);
						
						i_spectra[0]=_mm_add_ps(i_temp[0],i_spectra[0]);
						i_spectra[1]=_mm_add_ps(i_temp[1],i_spectra[1]);
						i_spectra[2]=_mm_add_ps(i_temp[2],i_spectra[2]);
						i_spectra[3]=_mm_add_ps(i_temp[3],i_spectra[3]);			
					} // for nTaps
					_mm_stream_ps(&spectra[(RpCl*c)*FpR + i*INNER*RpCl*FpR + bl*nChannels],i_spectra[0]);
					_mm_stream_ps(&spectra[(RpCl*c+1)*FpR + i*INNER*RpCl*FpR + bl*nChannels],i_spectra[1]);
					_mm_stream_ps(&spectra[(RpCl*c+2)*FpR + i*INNER*RpCl*FpR + bl*nChannels],i_spectra[2]);
					_mm_stream_ps(&spectra[(RpCl*c+3)*FpR + i*INNER*RpCl*FpR + bl*nChannels],i_spectra[3]);
				}//inner
			}// for nBlocks
		}// Outer
	}// parallel block
}




// Mk4 variant for AVX processors
void FIR_para_INT_Mk4_256(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nBlocks, int nThreads){
	int c,t,bl,th_id,block_step;
	nChannels=nChannels*2;
	
	//registers
	int FpR=8,RpCl=2;
	int REG=nChannels/(RpCl*FpR);
	__m256 i_data[2];// because of the size of the cacheline 64byte
	__m256 i_coeff[2];
	__m256 i_spectra[2];
	__m256 i_temp[2];
	
	omp_set_num_threads(nThreads);
	#pragma omp parallel shared(input_data,spectra,coeff) private(bl,t,c,i_data,i_coeff,i_spectra,i_temp)
	{
		//th_id = omp_get_thread_num();
		//nThreads = omp_get_num_threads();
		block_step=nBlocks/nThreads;
		
		#pragma omp for schedule(static,block_step)
		for(bl=0;bl<nBlocks-(nTaps-1);bl++){
			for(c=0;c<REG;c++){
				// Zeroing!
				i_spectra[0]=_mm256_setzero_ps();i_spectra[1]=_mm256_setzero_ps();
				
				for(t=0;t<nTaps;t++){
					//Loading data
					i_coeff[0]=_mm256_load_ps(&coeff[(RpCl*c)*FpR+t*nChannels]);
					i_coeff[1]=_mm256_load_ps(&coeff[(RpCl*c+1)*FpR+t*nChannels]);
										
					i_data[0]=_mm256_load_ps(&input_data[(t+bl)*nChannels+(RpCl*c)*FpR]);
					i_data[1]=_mm256_load_ps(&input_data[(t+bl)*nChannels+(RpCl*c+1)*FpR]);
					
					i_temp[0]=_mm256_mul_ps(i_data[0],i_coeff[0]);
					i_temp[1]=_mm256_mul_ps(i_data[1],i_coeff[1]);
					
					i_spectra[0]=_mm256_add_ps(i_temp[0],i_spectra[0]);
					i_spectra[1]=_mm256_add_ps(i_temp[1],i_spectra[1]);
				} // for nTaps
				_mm256_stream_ps(&spectra[(RpCl*c)*FpR + bl*nChannels],i_spectra[0]);
				_mm256_stream_ps(&spectra[(RpCl*c+1)*FpR + bl*nChannels],i_spectra[1]);
			}// nChannels
		}// for nBlocks
		
	}// parallel block	
}



void FIR_serial_INT_Mk4_256(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nBlocks){
	int c,t,bl,th_id,block_step;
	nChannels=nChannels*2;
	
	//registers
	int FpR=8,RpCl=2;
	int REG=nChannels/(RpCl*FpR);
	__m256 i_data[2];// because of the size of the cacheline 64byte
	__m256 i_coeff[2];
	__m256 i_spectra[2];
	__m256 i_temp[2];
	
		for(bl=0;bl<nBlocks-(nTaps-1);bl++){
			for(c=0;c<REG;c++){
				// Zeroing!
				i_spectra[0]=_mm256_setzero_ps();i_spectra[1]=_mm256_setzero_ps();
				
				for(t=0;t<nTaps;t++){
					//Loading data
					
					i_coeff[0]=_mm256_load_ps(&coeff[(RpCl*c)*FpR+t*nChannels]);
					i_coeff[1]=_mm256_load_ps(&coeff[(RpCl*c+1)*FpR+t*nChannels]);
										
					i_data[0]=_mm256_load_ps(&input_data[(t+bl)*nChannels+(RpCl*c)*FpR]);
					i_data[1]=_mm256_load_ps(&input_data[(t+bl)*nChannels+(RpCl*c+1)*FpR]);
					
					i_temp[0]=_mm256_mul_ps(i_data[0],i_coeff[0]);
					i_temp[1]=_mm256_mul_ps(i_data[1],i_coeff[1]);
					
					i_spectra[0]=_mm256_add_ps(i_temp[0],i_spectra[0]);
					i_spectra[1]=_mm256_add_ps(i_temp[1],i_spectra[1]);
				} // for nTaps
				_mm256_stream_ps(&spectra[(RpCl*c)*FpR + bl*nChannels],i_spectra[0]);
				_mm256_stream_ps(&spectra[(RpCl*c+1)*FpR + bl*nChannels],i_spectra[1]);
			}// nChannels
		}// for nBlocks
}


//**************
// PARALLEL CODE WITH INTRINSIC Mk4 BASED ON SERIAL Mk1 for architecture with register size = 128bit
//**************
//parallel version of the cache optimized code
void FIR_para_INT_Mk6_256(float *input_data, float *spectra, float *coeff, int nTaps, int nChannels, int nBlocks,int nThreads){
	//number of channels must by dividable by factor of 256
	int i,c,t,bl,th_id,block_step;
	int start,end;
	nChannels=nChannels*2;
	
	//registers
	int FpR=8,RpCl=2;
	int REG=nChannels/(RpCl*FpR);
	__m256 i_data[2];// because of the size of the cacheline 64byte
	__m256 i_coeff[2];
	__m256 i_spectra[2];
	__m256 i_temp[2];
	
	// Outer for
	// Inner for
	int INNER=16;
	int OUTER=REG/INNER;
	omp_set_num_threads(nThreads);
	#pragma omp parallel shared(input_data,spectra,coeff) private(start,end,i,bl,t,c,i_data,i_coeff,i_spectra,i_temp)
	{
		th_id = omp_get_thread_num();
		block_step=nBlocks/nThreads;
		if (block_step==0){
			if ( th_id<(nBlocks-(nTaps-1)) ) {
				start=th_id;
				end=th_id+1;
			}
		}
		else {
			end=(th_id+1)*block_step;
			if (th_id==nThreads-1) end=end-(nTaps-1);
			start=th_id*block_step;
		}
		
		for(i=0;i<OUTER;i++){
		
			for(bl=start;bl<end;bl++){
				for(c=0;c<INNER;c++){
					// Zeroing!
					i_spectra[0]=_mm256_setzero_ps();i_spectra[1]=_mm256_setzero_ps();
					
					for(t=0;t<nTaps;t++){
						//Loading data
						i_coeff[0]=_mm256_load_ps(&coeff[(RpCl*c)*FpR+t*nChannels]);
						i_coeff[1]=_mm256_load_ps(&coeff[(RpCl*c+1)*FpR+t*nChannels]);
											
						i_data[0]=_mm256_load_ps(&input_data[(t+bl)*nChannels+(RpCl*c)*FpR]);
						i_data[1]=_mm256_load_ps(&input_data[(t+bl)*nChannels+(RpCl*c+1)*FpR]);
						
						i_temp[0]=_mm256_mul_ps(i_data[0],i_coeff[0]);
						i_temp[1]=_mm256_mul_ps(i_data[1],i_coeff[1]);
						
						i_spectra[0]=_mm256_add_ps(i_temp[0],i_spectra[0]);
						i_spectra[1]=_mm256_add_ps(i_temp[1],i_spectra[1]);	
					} // for nTaps
					_mm256_stream_ps(&spectra[(RpCl*c)*FpR + bl*nChannels],i_spectra[0]);
					_mm256_stream_ps(&spectra[(RpCl*c+1)*FpR + bl*nChannels],i_spectra[1]);
				}//inner
			}// for nBlocks
		}// Outer
	}// parallel block
}
