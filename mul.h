#include "fftw3.h"
void mul(fftw_complex *fft_out_data,fftw_complex *ref){
	double temp[3200];
	for(int i=0;i<3200;i++){
		temp[i]=fft_out_data[i][0]*ref[i][0]-fft_out_data[i][1]*ref[i][1];
		fft_out_data[i][1]=ref[i][0]*fft_out_data[i][1]+fft_out_data[i][0]*ref[i][1];
		fft_out_data[i][0]=temp[i];
	}
}