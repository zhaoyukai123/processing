#include "stdlib.h"
void tran(double *a,int x0,int y0,int x1,int y1,double *temp){
	if((x1-x0+1)%2==0&&(x1-x0)>=(y1-y0)){
		tran(a,x0,y0,x1/2,y1,temp);
		tran(a,x1/2,y0,x1,y1,temp);
	}
	else if((y1-y0+1)%2==0&&(x1-x0)<=(y1-y0)){
		tran(a,x0,y0,x1,y1/2,temp);
		tran(a,x0,y1/2,x1,y1,temp);
	}
	else{

		for (int j=y0;j<y1+1;j++){
			for(int i=x0;i<x1+1;i++){
				temp[i*3200+j]=a[i+j*3200];
			}
		}
	}
	a=temp;
}
