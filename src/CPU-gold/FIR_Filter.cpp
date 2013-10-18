//TODO! if I will keep the current buffer arangement with blocks of size nTaps I need to change coeff to be more compatible with the chosen scheme
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
#include <math.h>
#include <fftw3.h>

#include "Filter_window.h"

using namespace std;


long int File_size_byte(ifstream &FILEIN) {
	ifstream::pos_type size;
	FILEIN.seekg(0,ios::end);
	size=FILEIN.tellg();
	FILEIN.seekg(0,ios::beg);
	return(size);
}

long int File_size_row(ifstream &FILEIN){
	std::size_t count=0;
	FILEIN.seekg(0,ios::beg);
	for(std::string line; std::getline(FILEIN, line); ++count){}
	return(count);
}


bool Load_data(ifstream &FILEIN,float *time, float *real, float *img, long int size){
	long int count=0;
	while (!FILEIN.eof()){
		FILEIN >> time[count] >> real[count] >> img[count];
		count++;
	}
	return(1);
}

void Create_buffer(float *time, float *real, float *img, float *wbuffer, unsigned int nTaps, unsigned int nSamples, unsigned int *oldestblock){
	//TODO! whole function must be changed when parallelized
	// this function is to be called only when filling the wbuufer for the first time
	// copying real part
	for(int s=0;s<nSamples;s++){
		for(int t=0;t<nTaps;t++){
			wbuffer[s*nTaps+t]=real[nSamples*t+s]
			wbuffer[s*nTaps+t+nTaps*nSamples]=img[nSamples*t+s];
		}
	}

	// TODO! need to be changed when parallelized
	*oldestblock=0;
}

void Update_buffer(float *time, float *real, float *img, float *wbuffer, unsigned int nTaps, unsigned int nSamples, long int blocktocopy,unsigned int *oldestblock){
	//memcpy(&wbuffer[oldestblock*nSamples],&real[blocktocopy*nSamples],nSamples);
	//memcpy(&wbuffer[oldestblock*nSamples+nSamples*nTaps],&img[blocktocopy*nSamples],nSamples);
	//TODO! change for parallelization
	for(int s=0;s<nSamples;s++){
		wbuffer[s*nTaps+*oldestsample]=real[blocktocopy*nSamples+s];
		wbuffer[s*nTaps+*oldestsample+nTaps*nSamples]=img[blocktocopy*nSamples+s];
	}
	
	*oldestblock++;
}

void load_buffer(){
	

} 

int main(void) {
	// definition of constants
	const int nTaps=8;
	const int nSamples=512;
	const int register_size=256;//CPU
	const int in_buffer_size=register_size/32;	
	const int w_buffer_size=nSamples*nTaps*3;

	float *time,*real,*img;
	float *coeff;
	long int filesize=0;
	bool flag=0;
	coeff= new float [nTaps*nSamples];
	
	// Loading data from file
	ifstream FILEIN;
	FILEIN.open("test.dat",ios::in);
	if (!FILEIN.fail()){
		filesize=File_size_row(FILEIN);
		time=new float[filesize];
		real=new float[filesize];
		img=new float[filesize];		
		FILEIN.clear();
		FILEIN.seekg(0,ios::beg);
		flag=Load_data(FILEIN,time,real,img,filesize);
		if (flag) { 
			cout << "File opened and loaded" << endl;
		}
	}
	else {
		cout << "File not found" << endl;
		exit(2);
	}
	FILEIN.close();
	
	if(filesize%nSamples!=0){
		printf("Error in input file! Data are not in multiple of nSamples!");
		exit(12);
	}

	unsigned int nBlocks=filesize/nSamples;

	//****************************************************************************
	// Testing if loading was succesful
	cout << "Number of rows: " << filesize << endl;
	for (int f=0;f<filesize;f++){
		cout << time[f] << " " << real[f] << " " << img[f] << endl;
	}	
	//****************************************************************************
	
	
	
	// geting window coefficients
	coeff = new float [nSamples * nTaps];
	//Creating filter window
	Window_coefficients(nTaps,nSamples,1,coeff);
	
	//****************************************************************************
	// Testing if there is something 
	for(int f=0;f<20;f++){
		cout << coeff[f] << " - ";
	}
	cout << endl;
	//****************************************************************************

	// Memory arrangement

	// I assume that format of the input signal is
	// | <----nSamples----> |              |           |          |   |             |
	// |____________________|______________|___________|__________|...|_____________|
	//    Oldest time block                                             newest time block
	
	
	long int Block=0;
	unsigned int Sample=0;
	unsigned int oldestblock=0;
	float *multi1_buffer,*multi2_buffer,*multi1_c_buffer,*multi2_c_buffer;
	float *add1_1_buffer,*add1_2_buffer,*add2_1_buffer,*add2_2_buffer;
	float *wbuffer;
	float *result;
	// creating buffers with insintric and parallelization in mind
	// CPU insintric: register size 256bit
	// float size: 32bit => buffer size is 8 floats
	// WARNING! in_buffer_size must be = to nTaps (at this stage)
	multi1_buffer=new float [in_buffer_size];
	multi2_buffer=new float [in_buffer_size];
	multi1_c_buffer=new float [in_buffer_size];
	multi2_c_buffer=new float [in_buffer_size];
	add1_1_buffer=new float [in_buffer_size];
	add1_2_buffer=new float [in_buffer_size];
	add2_1_buffer=new float [in_buffer_size];
	add2_2_buffer=new float [in_buffer_size];
	
	wbuffer=new float [w_buffer_size];

	result=new float [nSamples];
	// implementation of this will differ depending on parallelization (?)
	// buffer structure
	// |___| <= one time block
	//  V   
	// | # # # # # # # | <= #=One time block number of # is nTaps
	//   V
	// |--*1--|--*2--|--*3--| <= x nSamples (*1 time *2 real *3 img)  for now without time!
	Create_buffer(time,real,img,wbuffer,nTaps,nSamples,&oldestblock);
	for(Block=nTaps;Block<nBlocks;Block++){
		for(Sample=0,Sample<nSample/2;Sample++){
			//WARNING these must be tail which will take care of the remaining rows in the insintric buffers 
			
			//loading data into multi buffers
			memcpy(multi1_buffer,&wbuffer[Sample*nTaps],nTaps);
			memcpy(multi2_buffer,&wbuffer[Sample*nTaps],nTaps);
			//loading coefficients into multi buffers
			//not done
			//not done should be through memcpy
			
			// this will in future be taken care of by insintric operation on registers (hopefully) 
			// this multiplies two multi insintric buffers doing c*x operation and copies the result into first half copydestination=0 or second half copydestination =1 add1_buffer for first stage addition
			multiply_by_coeff(multi1_buffer,multi1_c_buffer,0,add1_1_buffer,add1_2_buffer);
			multiply_by_coeff(multi2_buffer,multi2_c_buffer,1,add1_1_buffer,add1_2_buffer);
			
			// this should be taken care of by insintric operation
			addition_stage1(add1_1_buffer,add1_2_buffer,add2_1_buffer,add2_2_buffer);
			
			addition_stage2(add2_1_buffer,add2_2_buffer);
			addition_stage3(add2_1_buffer,add2_2_buffer,result,Sample);
			copy_result();//this copies the result to the final array representing one spectra
		}
		//TODO! Add tail
		Update_buffer(time,real,img,wbuffer,nTaps,nSamples,Block,&oldestblock);
	}




	
	delete[] time;
	delete[] real;
	delete[] img;

	delete[] coeff;

	delete[] multi1_buffer;
	delete[] multi2_buffer;
	delete[] add1_buffer;
	delete[] add2_buffer;

    #pragma omp parallel
    printf("Hello, world.\n");
    system("PAUSE");
    return EXIT_SUCCESS;
}

