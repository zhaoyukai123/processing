#include<math.h>
#include<omp.h>
#include <complex>
using namespace std;
void FFT1(complex<double> *raw_data_fft,int n,int k,complex<double> *raw_data_fft_result,int l,int il) 
{	
	int it,m,is,i,j,nv,l0;
	double p,q,s,vr,vi,poddr,poddi;
	for (it=0; it<=n-1; it++)
	{ 
		m=it; is=0;
		for (i=0; i<=k-1; i++)
		{ 
			j=m/2; 
			is=2*is+(m-2*j);
			m=j;
		}
		raw_data_fft_result[it].real(raw_data_fft[is].real());
		raw_data_fft_result[it].imag(raw_data_fft[is].imag());
	}
	raw_data_fft[0].real(1.0); 
	raw_data_fft[0].imag(0.0);
	p=6.283185306/(1.0*n);
	raw_data_fft[1].real(cos(p));
	raw_data_fft[1].imag(-sin(p));
	if (l!=0) 
		raw_data_fft[1].imag(-raw_data_fft[1].imag());
	for (i=2; i<=n-1; i++)
	{ 
		p=raw_data_fft[i-1].real()*raw_data_fft[1].real();
		q=raw_data_fft[i-1].imag()*raw_data_fft[1].imag();
		s=(raw_data_fft[i-1].real()+raw_data_fft[i-1].imag())*(raw_data_fft[1].real()+raw_data_fft[1].imag());
		raw_data_fft[i].real(p-q); 
		raw_data_fft[i].imag(s-p-q);
	}
	for (it=0; it<=n-2; it=it+2)
	{ 
		vr=raw_data_fft_result[it].real(); 
		vi=raw_data_fft_result[it].imag();
		raw_data_fft_result[it].real(vr+raw_data_fft_result[it+1].real()); 
		raw_data_fft_result[it].imag(vi+raw_data_fft_result[it+1].imag());
		raw_data_fft_result[it+1].real(vr-raw_data_fft_result[it+1].real()); 
		raw_data_fft_result[it+1].imag(vi-raw_data_fft_result[it+1].imag());
	}
	m=n/2;
	nv=2;
	for (l0=k-2; l0>=0; l0--)
	{ 
		m=m/2; 
		nv=2*nv;
		for (it=0; it<=(m-1)*nv; it=it+nv)
			for (j=0; j<=(nv/2)-1; j++)
			{ 
				p=raw_data_fft[m*j].real()*raw_data_fft_result[it+j+nv/2].real();
				q=raw_data_fft[m*j].imag()*raw_data_fft_result[it+j+nv/2].imag();
				s=raw_data_fft[m*j].real()+raw_data_fft[m*j].imag();
				s=s*(raw_data_fft_result[it+j+nv/2].real()+raw_data_fft_result[it+j+nv/2].imag());
				poddr=p-q;
				poddi=s-p-q;
				raw_data_fft_result[it+j+nv/2].real(raw_data_fft_result[it+j].real()-poddr);
				raw_data_fft_result[it+j+nv/2].imag(raw_data_fft_result[it+j].imag()-poddi);
				raw_data_fft_result[it+j].real(raw_data_fft_result[it+j].real()+poddr);
				raw_data_fft_result[it+j].imag(raw_data_fft_result[it+j].imag()+poddi);
			}
	}
	if (l!=0)
		for (i=0; i<=n-1; i++)
		{ 
			raw_data_fft_result[i].real(raw_data_fft_result[i].real()/(1.0*n));
			raw_data_fft_result[i].imag(raw_data_fft_result[i].imag()/(1.0*n));
		}
	if (il!=0)
		for (i=0; i<=n-1; i++)
		{ 
			raw_data_fft[i].real(sqrt(raw_data_fft_result[i].real()*raw_data_fft_result[i].real()+raw_data_fft_result[i].imag()*raw_data_fft_result[i].imag()));
			if (fabs(raw_data_fft_result[i].real())<0.000001*fabs(raw_data_fft_result[i].imag()))
			{ 
				if ((raw_data_fft_result[i].imag()*raw_data_fft_result[i].real())>0)
						raw_data_fft[i].imag(90.0);
				else raw_data_fft[i].imag(-90.0);
			}
			else
				raw_data_fft[i].imag(atan(raw_data_fft_result[i].imag()/raw_data_fft_result[i].real())*360.0/6.283185306);
		}

	return;
}