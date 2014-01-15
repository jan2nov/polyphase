#include <iostream>
#include <fstream>
#include <iomanip> 

using namespace std;

long int File_size_row(ifstream &FILEIN){
		std::size_t count=0;
		FILEIN.seekg(0,ios::beg);
		for(std::string line; std::getline(FILEIN, line); ++count){}
	return((long int)count);
}

bool Load_data(ifstream &FILEIN, float2 *data){
	long int count = 0;
	while (!FILEIN.eof()){
		FILEIN >> data[count].x;
		data[count].y = 0;
		count++;
	}
	return(1);
}

bool Save_data_spectra(char str[],float2 *spectra, unsigned int size){
	ofstream FILEOUT;
	FILEOUT.open(str);
	for (unsigned int c = 0; c < size; c++){
		FILEOUT << spectra[c].x << " " << spectra[c].y << endl;
	}
	FILEOUT.close();
	return(1);
}


bool Save_window_data(float *coeff, unsigned int nChannels, unsigned int nTaps){
	ofstream FILEOUT;
	FILEOUT.open("window.tap");

	for(unsigned int t=0;t<nTaps;t++){
		for (unsigned int c=0;c<nChannels;c++){
			FILEOUT << coeff[t*nChannels+c] << endl;
		}
	}

	FILEOUT.close();
	return (true);
}

bool Load_window_data(float *coeff){
	unsigned int count=0;
	ifstream FILEIN;
	FILEIN.open("window.tap",ios::in);
	while (!FILEIN.eof()){
		FILEIN >> coeff[count];
		count++;
	}
	return(1);
}

bool save_time(char str[], int num_blocks, float fir_time, float fft_time, float mem_time_in, float mem_time_out, const int nChannels, const int nTaps){
	ofstream FILEOUT;
	FILEOUT.open (str, std::ofstream::out | std::ofstream::app);
	double flops = 0.0f; 
	double bandwidth = 0.0f;
	flops = (3.0*nChannels*nTaps*sizeof(float)+nChannels*2*sizeof(float))*(num_blocks)*1000/fir_time;
	bandwidth = (4.0*nTaps)*nChannels*(num_blocks)*1000.0/fir_time;
	//------------------
		FILEOUT << std::fixed << std::setprecision(8) << num_blocks << "\t" << fir_time/1000 << "\t" << (fir_time+fft_time)/1000 << "\t" << std::scientific << flops << "\t" << bandwidth << "\t" << std::fixed << mem_time_in/1000  << "\t" << mem_time_out/1000  << endl;
	FILEOUT.close();
	return 0;
}

