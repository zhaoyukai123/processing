#include "stdio.h"
#include "stdlib.h"
#include <math.h>

void circshift(double *cor_data,int *shift,int length){
	double temp[6400];
	for (int i=0;i<length;i++){
		for(int j=0;j<3200;j++){
			temp[j+3200]=cor_data[j+i*3200];
			temp[j]=cor_data[j+i*3200];
		}
		for(int j=0;j<3200;j++){
			cor_data[j+i*3200]=temp[j+shift[i]];
		}
	}
}
