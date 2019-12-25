#include "stdlib.h"
void cor(double cor_data[][1000],double in[][10000]){
	double temp[10000][3200];
	for (int i=0;i<10000;i++){
		for (int j=0;j<3200;j++){
			temp[i][j]=in[j][i];
	}
	for(int i=0;i<1000;i++){
		for(int j=0;j<3200;j++){
			temp[i*10][j]=(temp[i*10][j]+temp[i*10+1][j]+temp[i*10+2][j]+temp[i*10+3][j]+temp[i*10+4][j]+temp[i*10+5][j]+temp[i*10+6][j]+temp[i*10+7][j]+temp[i*10+8][j]+temp[i*10+9][j])/10;
		}
	}
	for (int i=0;i<1000;i++){
		for (int j=0;j<3200;j++){
			cor_data[j][i]=temp[i][j];
		}
	}
}
}