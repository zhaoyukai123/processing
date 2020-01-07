#include "stdio.h"
#include "stdlib.h"
#include <math.h>
#include "fliplr.h"
#include "BLEndianUshort.h"
#include "circshift.h"
#include "time.h"
#include <omp.h>
#include "cor.h"
#include "incor.h"
#include "lslf.h"
#define M_PI 3.14159265358979323846
typedef struct complex{
    double real;
    double imag;
};
complex add(complex num1, complex num2);
complex subtract(complex num1, complex num2);
complex multiply(complex num1, complex num2);
complex divide(complex num1, complex num2);
complex power(complex num, double n);
complex complex_exp(double theta);
double magnitude(complex num);
double phase(complex num);
void fft(complex input[], complex * output[], int n);
complex conjugate(complex num);
void fft_driver(complex input[], complex output[], int n, int step){
    complex diff, cexp;

    if (step < n) {

        fft_driver(output, input, n, step * 2);
        fft_driver(output + step, input + step, n, step * 2);

        for (int i = 0; i < n; i += 2 * step) {
            diff = multiply(complex_exp(-M_PI * (double) i / (double) n), output[i + step]);
            input[i / 2] = add(output[i], diff);
            input[(i + n) / 2] = subtract(output[i], diff);
        }
    }
}

void fft(complex *in,complex *out,complex *dummy,int n){
	for (int i = 0; i < n; ++i) {
    dummy[i]=in[i];
    out[i]=in[i];
	}
	fft_driver(out, dummy, n, 1);
}
void ifft_driver(complex input[], complex output[], int n, int step){
    complex diff, cexp;

    if (step < n) {

        fft_driver(output, input, n, step * 2);
        fft_driver(output + step, input + step, n, step * 2);

        for (int i = 0; i < n; i += 2 * step) {
            diff = multiply(complex_exp(M_PI * (double) i / (double) n/3200), output[i + step]);
            input[i / 2] = add(output[i], diff);
            input[(i + n) / 2] = subtract(output[i], diff);
        }
    }
}

void ifft(complex *in,complex *out,complex *dummy,int n){
	for (int i = 0; i < n; ++i) {
    dummy[i]=in[i];
    out[i]=in[i];
	}
	ifft_driver(out, dummy, n, 1);
}
void main(){
	double *a = NULL;
	int n=1000000000;
	int top,bottom,middle;
	int error=0;
	while (error==0){
    a = (double*)malloc( n * sizeof(double)*13076);
    if (!a){
		free(a);
		top=n;
		n/=2;
		bottom=n;
    }
	else error=1;
}
	while(error==1){
		middle=top/2+bottom/2;
		a = (double*)malloc( middle* sizeof(double)*13076);
    if (!a){
		free(a);
		top=middle;
		}
	else if(top-bottom>500){
		bottom=middle;	
		free(a);
		}
	else{
		free(a);
		n=middle;
		error=0;
	}

}
	n=n/100000;
	n=n-n%50;
	if(n>traces)n=traces-traces%50;
	int data_col=3200;
	FILE *bxds;
	bxds = fopen("d:/position/bxds","rb");
	FILE *bxds1;
	bxds1 = fopen("d:/position/bxds1","wb");
	fseek(bxds,0L,SEEK_END);
	long ENum=ftell(bxds);
	fseek(bxds,0L,SEEK_SET);
	long SNum=ftell(bxds);
	long ByteNum=ENum-SNum;
	printf("%lu\n",ByteNum);
	int traces;
	traces=ByteNum/12870;
	printf("%d\n",traces);
	double I[100]={-63, -92, -109, -75, -87, -50, -116, -154, -22, 68, 141,-610, 1461, 3807, -6147, -5375, 10651, -4412, -9810, 15386,-3070, -14499, 15130, 3677, -15935, 3743, 13362, -5884, -13301, 8455, 12542, -8744, -11977, 5105, 13754, -961,-14342, -5184, 10294, 12194, -2709, -14352, -8807, 5965,15350, 8368, -6605, -14990, -11515, 196, 11276, 15490,10300, -645, -10730, -15307, -13379, -7342, 377, 7264,11662, 13435, 11530, 3243, -4865, -4427, -3233, -4000,-2472, -2498, -2361, -1230, -1311, -618, -578, -569, -121,-319, 206, 328, 436, 613, 318, 514, 353, 277, 221, 34, 250,132, 199, 189, 75, 190, 65, 106, 19, -64, -14, -117};
	fliplr(I,100);
	complex *in,*out, *dummy;
	in  = (complex*)malloc(sizeof(complex) * 3200);
	for(int i=0;i<3200;i++){
		in[i].real=0;	
		in[i].imag=0;
	}
	for(int i=0;i<100;i++){
		in[i].real=I[i];	
	}
	double x[3200];
	for(int i=0;i<3200;i++)x[i]=i;
	out = (complex*)malloc(sizeof(complex) * 3200);
	fft(in,out,dummy,3200);
	int max_freq=1120;
	int min_freq=160;
	int dfreq=961;
	double linspace[961];
	double hamming[961];
	for(int i=0;i<961;i++){
		linspace[i]=i*M_PI/961;
		hamming[i]=sin(linspace[i]);
	}
	double hfilter[3200];
	for(int i=0;i<3200;i++){
		hfilter[i]=0;
	}
	for(int i=160;i<1121;i++){
		hfilter[i]=hamming[i-160]*2;
	}
	for(int i=0;i<3200;i++){
		out[i].real=out[i].real*hfilter[i];
		out[i].imag=out[i].imag*hfilter[i];
	}
	//process
	complex *fft_in_data,*fft_out_data;
	fft_in_data  = (complex*) malloc(sizeof(complex) * 3200);
	fft_out_data  = (complex*) malloc(sizeof(complex) * 3200);
	double *in_data,*cor_data,*incor_data,*result;
	in_data=(double*) malloc(sizeof(double) * 3200*n);
	cor_data=(double*) malloc(sizeof(double) * 3200*n/10);
	incor_data=(double*) malloc(sizeof(double) * 3200*n/50);
	result=(double*) malloc(sizeof(double) * 3200*n/50);
	short *temp;
	temp=(short*) malloc(sizeof(short) * 3200);
	int *max_temp,*shift;
	max_temp=(int*) malloc(sizeof(int) * n/10);
	max_temp=(int*) malloc(sizeof(int) * n/10);
	int read_traces=n;
	int read_col=3200;
	#pragma omp parallel for
	for(int b=0;b<traces/n;b++){


	for(int i=0;i<n;i++){
		fseek(bxds,70*(i+1)+i*12800+b*12870*10000,0);
		fread(temp,sizeof(short),3200,bxds);
		for(int j=0;j<3200;j++){
			temp[j]=BLEndianUshort(temp[j]);
			in_data[i*3200+j]=temp[j];
		}
	}
	double *temp1,*temp2,*temp3;
	temp1=(double*)malloc(sizeof(double)*n);
	temp2=(double*)malloc(sizeof(double)*n/10);
	temp3=(double*)malloc(sizeof(double)*n/50);
	cor(cor_data,in_data,temp1,temp2,n);
	//for(int i=0;i<read_traces/10;i++){
	//	for(int j=0;j<3200;j++){
	//		cor_data[j][i]=(in_data[j][i*10]+in_data[j][i*10+1]+in_data[j][i*10+2]+in_data[j][i*10+3]+in_data[j][i*10+4]+in_data[j][i*10+5]+in_data[j][i*10+6]+in_data[j][i*10+7]+in_data[j][i*10+8]+in_data[j][i*10+9])/10;
	//	}
	//}

	for(int i=0;i<read_traces/10;i++){
		max_temp[i]=-99999;
	}
	for(int i=0;i<n/10;i++){
		for(int j=0;j<3200;j++){
			if(cor_data[j+i*3200]>max_temp[i]){
				max_temp[i]=cor_data[j+i*3200];
				shift[i]=j;
			}
		}
	}

	circshift(cor_data,shift,n/10);
	lslf(cor_data,x);

	for (int i=0;i<n/10;i++){
		for(int j=0;j<3200;j++){
			fft_in_data[j].real=cor_data[j+i*3200];
			fft_in_data[j].imag=0;
		}
		fft(fft_in_data,fft_out_data,dummy,3200);
		for(int j=0;j<3200;j++){
			multiply(fft_out_data[j],out[j]);
			conjugate(fft_out_data[j]);
	}	
		ifft(fft_out_data,fft_in_data,dummy,3200);

		for(int j=0;j<3200;j++){
			cor_data[j+i*3200]=magnitude(fft_in_data[j]);;
		}
	}
	for(int i=0;i<n/10;i++){
		shift[i]=3200-shift[i];
	}
	circshift(cor_data,shift,n/10);
	incor(incor_data,cor_data,temp3,temp3,n);
	//	for(int i=0;i<read_traces/50;i++){
	//		for(int j=0;j<3200;j++){
	//			incor_data[j][i]=(cor_data[j][i*5]/5+cor_data[j][i*5+1]/5+cor_data[j][i*5+2]/5+cor_data[j][i*5+3]/5+cor_data[j][i*5+4]/5);
	//	}
	//}
		for(int i=0;i<n/50;i++){
			for(int j=0;j<3200;j++){
				incor_data[j+i*3200]=20000*log10(incor_data[j+i*3200]);
		}
	}
	
	for(int i=0;i<n/50;i++){
			for(int j=0;j<3200;j++){
				result[i*3200+j]=incor_data[j+i*3200];
		}
	}

	fwrite(result,sizeof(double),3200*n/50,bxds1);
}
}
double magnitude(complex num) {
    return sqrt(pow(num.real, 2) + pow(num.imag, 2));
}

double phase(complex num) {
    return atan(num.imag / num.real);
}

complex add(complex num1, complex num2) {
    complex result;
    result.real = num1.real + num2.real;
    result.imag = num1.imag + num2.imag;

    return result;
}

complex conjugate(complex num) {
    complex result;
    result.real = num.real;
    result.imag = -num.imag;

    return result;
}

complex subtract(complex num1, complex num2) {
    complex result;
    result.real = num1.real - num2.real;
    result.imag = num1.imag - num2.imag;

    return result;
}

complex multiply(complex num1, complex num2) {
    complex result;
    result.real = num1.real * num2.real - num1.imag * num2.imag;
    result.imag = num1.real * num2.imag + num1.imag * num2.real;

    return result;
}

complex divide(complex num1, complex num2) {
    complex result;

    double magnitude_square = pow(num2.real, 2) + pow(num2.imag, 2);
    num2.imag = -num2.imag;

    result = multiply(num1, num2);
    result.real = result.real / magnitude_square;
    result.imag = result.imag / magnitude_square;

    return result;
}

complex power(complex num, double n) {
    complex result;

    double magnitude = sqrt(pow(num.real, 2) + pow(num.imag, 2));

    result.real = pow(magnitude, n) * cos(n * atan(num.imag / num.real));
    result.imag = pow(magnitude, n) * sin(n * atan(num.imag / num.real));

    return result;
}
complex complex_exp(double theta) {
    complex result;
    result.real = cos(theta);
    result.imag = sin(theta);

    return result;
}
