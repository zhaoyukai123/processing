#include "fftw3.h"
#include "math.h"
void abs(fftw_complex *in,fftw_complex *out){
	for(int i=0;i<3200;i++){
		out[i][0]=sqrt(in[i][0]*in[i][0]+in[i][1]*in[i][1]);
	}
}