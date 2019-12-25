#include "stdlib.h"
void incor(double cor_data[][200],double in[][1000]){
	double temp[1000][3200];
	for (int i=0;i<1000;i++){
		for (int j=0;j<3200;j++){
			temp[i][j]=in[j][i];
	}
	for(int i=0;i<2000;i++){
		for(int j=0;j<3200;j++){
			temp[i*10][j]=(temp[i*10][j]+temp[i*10+1][j]+temp[i*10+2][j]+temp[i*10+3][j]+temp[i*10+4][j])/5;
		}
	}
	for (int i=0;i<200;i++){
		for (int j=0;j<3200;j++){
			cor_data[j][i]=temp[i][j];
		}
	}
}
}