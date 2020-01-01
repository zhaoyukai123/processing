#include "stdlib.h"
#include "tran.h"
void cor(double *cor_data,double *in,double *temp1,double *temp2,int n){
	tran(in,0,0,n,3200,temp1);
	for(int i=0;i<n/10;i++){
		for(int j=0;j<3200;j++){
			cor_data[i*3200+j]=(temp1[i*10+j]+temp1[i*10+1+j]+temp1[i*10+2+j]+temp1[i*10+3+j]+temp1[i*10+4+j]+temp1[i*10+5+j]+temp1[i*10+6+j]+temp1[i*10+7+j]+temp1[i*10+8+j]+temp1[i*10+9+j])/10;
		}
	}
	tran(cor_data,0,0,n/10,3200,temp2);
}
