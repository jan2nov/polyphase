#include <iostream>
#include <fstream>

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
