#include "stdlib.h"
void incor(double *cor_data,double *in,double *temp2,double *temp3,int n){
	tran(in,0,0,n/10,3200,temp2);
	for(int i=0;i<n/50;i++){
		for(int j=0;j<3200;j++){
			cor_data[i*3200+j]=(temp2[i*5+j]+temp2[i*5+1+j]+temp2[i*5+2+j]+temp2[i*5+3+j]+temp2[i*5+4+j])/5;
		}
	}
	tran(cor_data,0,0,n/50,3200,temp3);
}
