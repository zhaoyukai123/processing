#include "stdlib.h"
void lslf(double *y,double *x)
{
	double a,b;
	double sum_x2,sum_y,sum_x,sum_xy;
for(int j=0;j<1000;j++){
	sum_x2 = 0.0;sum_y  = 0.0;sum_x  = 0.0;sum_xy = 0.0;
	for (int i = 0; i < 3200; ++i) {
	sum_x2 += x[i]*x[i];
	sum_y+=y[i+j*3200];
	sum_x+=x[i];
	sum_xy += x[i]*y[i+j*3200];
}
	a = (3200*sum_xy - sum_x*sum_y)/(3200*sum_x2 - sum_x*sum_x);
    b = (sum_x2*sum_y - sum_x*sum_xy)/(3200*sum_x2-sum_x*sum_x);
	for (int i = 0; i < 3200; ++i) {
		y[i+j*3200]-=(a*x[i]+b);
	}
}
}
