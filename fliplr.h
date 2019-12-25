#include "stdio.h"
#include "stdlib.h"
#include <math.h>
#define M_PI 3.14159265358979323846

void fliplr(double *I,int length){
	for(int i=0;i<length/2;i++)
	{
		double temp;
		temp=I[i];
		I[i]=I[length-1-i];
		I[length-1-i]=temp;
	}
}