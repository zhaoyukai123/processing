#include "stdio.h"
#include "stdlib.h"
#include <math.h>
#include "fliplr.h"
#include "fftw3.h"
#include "BLEndianUshort.h"
#include "circshift.h"
#include "mul.h"
#include "abs.h"
#include "time.h"
#include <omp.h>
#include "lslf.h"
#include "cor.h"
#include "incor.h"
#define M_PI 3.14159265358979323846

void main(){
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
	fftw_complex *out,*in;
	in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 3200);
	for(int i=0;i<3200;i++){
		in[i][0]=0;	
		in[i][1]=0;
	}
	for(int i=0;i<100;i++){
		in[i][0]=I[i];	
	}
	double x[3200];
	for(int i=0;i<3200;i++)x[i]=i;
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 3200);
	fftw_plan p = fftw_plan_dft_1d(3200, in, out, FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	//fftw_free(out);
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
		out[i][0]=out[i][0]*hfilter[i];
		out[i][1]=out[i][1]*hfilter[i];
	}
	//process
	fftw_complex *fft_in_data,*fft_out_data;
	fft_in_data  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 3200);
	fft_out_data  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 3200);
	double in_data[3200][10000];
	double cor_data[3200][10000/10];
	short temp[3200];
	int max_temp[10000/10];
	int shift[10000/10];

	int read_traces=10000;
	int read_col=3200;
	double incor_data[3200][10000/50];
	double result[10000/50][3200];
	#pragma omp parallel for
	for(int b=0;b<traces/10000;b++){


	for(int i=0;i<read_traces;i++){
		fseek(bxds,70*(i+1)+i*12800+b*12870*10000,0);
		fread(temp,sizeof(short),3200,bxds);
		for(int j=0;j<3200;j++){
			temp[j]=BLEndianUshort(temp[j]);
			in_data[j][i]=temp[j];
		}
	}
	cor(cor_data,in_data);
	//for(int i=0;i<read_traces/10;i++){
	//	for(int j=0;j<3200;j++){
	//		cor_data[j][i]=(in_data[j][i*10]+in_data[j][i*10+1]+in_data[j][i*10+2]+in_data[j][i*10+3]+in_data[j][i*10+4]+in_data[j][i*10+5]+in_data[j][i*10+6]+in_data[j][i*10+7]+in_data[j][i*10+8]+in_data[j][i*10+9])/10;
	//	}
	//}

	for(int i=0;i<read_traces/10;i++){
		max_temp[i]=-99999;
	}
	for(int i=0;i<read_traces/10;i++){
		for(int j=0;j<3200;j++){
			if(cor_data[j][i]>max_temp[i]){
				max_temp[i]=cor_data[j][i];
				shift[i]=j;
			}
		}
	}

	circshift(cor_data,shift,read_traces/10);
	lslf(cor_data,x);

	for (int i=0;i<read_traces/10;i++){
		for(int j=0;j<3200;j++){
			fft_in_data[j][0]=cor_data[j][i];
			fft_in_data[j][1]=0;
		}
		fftw_plan q,r;
		q = fftw_plan_dft_1d(3200, fft_in_data, fft_out_data, FFTW_FORWARD,FFTW_ESTIMATE);
		fftw_execute(q);
		fftw_destroy_plan(q);
		mul(fft_out_data,out);
		r = fftw_plan_dft_1d(3200, fft_out_data, fft_in_data, FFTW_BACKWARD,FFTW_ESTIMATE);
		fftw_execute(r);
		fftw_destroy_plan(r);
		abs(fft_in_data,fft_out_data);
		for(int j=0;j<3200;j++){
			cor_data[j][i]=fft_out_data[j][0];
		}
	}
	for(int i=0;i<read_traces/10;i++){
		shift[i]=3200-shift[i];
	}
	circshift(cor_data,shift,read_traces/10);
	incor(incor_data,cor_data);
	//	for(int i=0;i<read_traces/50;i++){
	//		for(int j=0;j<3200;j++){
	//			incor_data[j][i]=(cor_data[j][i*5]/5+cor_data[j][i*5+1]/5+cor_data[j][i*5+2]/5+cor_data[j][i*5+3]/5+cor_data[j][i*5+4]/5);
	//	}
	//}
		for(int i=0;i<read_traces/50;i++){
			for(int j=0;j<3200;j++){
				incor_data[j][i]=20000*log10(incor_data[j][i]);
		}
	}
	
	for(int i=0;i<read_traces/50;i++){
			for(int j=0;j<3200;j++){
				result[i][j]=incor_data[j][i];
		}
	}

	fwrite(result,sizeof(double),3200*read_traces/50,bxds1);
}
}
