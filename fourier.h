#include <math.h>	// sqrt

#define PI 3.14159265359
#define TWOPI	(2.0*PI)

void four1(double data[], int nn, int isign)
{
    int n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
    
    n = nn << 1;
    j = 1;
    for (int i = 1; i < n; i += 2) {
	if (j > i) {
	    tempr = data[j];     data[j] = data[i];     data[i] = tempr;
	    tempr = data[j+1]; data[j+1] = data[i+1]; data[i+1] = tempr;
	}
	m = n >> 1;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }
    mmax = 2;
    while (n > mmax) {
	istep = 2*mmax;
	theta = TWOPI/(isign*mmax);
	wtemp = sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j =i + mmax;
		tempr = wr*data[j]   - wi*data[j+1];
		tempi = wr*data[j+1] + wi*data[j];
		data[j]   = data[i]   - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr)*wpr - wi*wpi + wr;
	    wi = wi*wpr + wtemp*wpi + wi;
	}
	mmax = istep;
    }
}

//zero padding of complex coefficients
void realft(double* data, unsigned long n)
{
	double *X = new double[2*n+1]; 
	for(int i=0; i<n; i++)
	{
		X[2*i+1] = data[i];
		X[2*i+2] = 0.0;
	}	
	
	four1(X, n, 1);
	
	//double *res = new double[n];
	
	for(int i=0; i<n; i++) data[i] = X[i+1];
}

//Calculates the Fourier transform of a set of n real-valued data points. 
//Replaces this data (which is stored in array data[1..n]) 
//by the positive frequency half of its complex Fourier transform.
//The real-valued first and last components of the complex transform are returned as elements
//data[1] and data[2], respectively. n must be a power of 2. This routine also calculates the
//inverse transform of a complex data array if it is the transform of real data. 
//(Result in this case must be multiplied by 2/n.)
/*
void realft(double* data, unsigned long n)
{
	unsigned long i,i1,i2,i3,i4,np3;
	double c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;
	//Double precision for the trigonometric recurrences.
	
	theta=3.141592653589793/(double) (n>>1);
	
	//Initialize the recurrence.
	c2 = -0.5;
	four1(data,n>>1);
	//The forward transform is here.

	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;

	for (i=1;i<(n>>2);i++) {
		//Case i=0 done separately below.
		//i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		i2=1+(i1=i+i);
		i4=1+(i3=n-i1);
		
		h1r=c1*(data[i1]+data[i3]);
		//The two separate transforms are separated out of data.
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		//Here they are recombined to form the true transform of the original real data.
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		//The recurrence.
		wi=wi*wpr+wtemp*wpi+wi;
	}

	data[0] = (h1r=data[0])+data[1];
	//Squeeze the first and last data together to get them all within the original array.
	data[1] = h1r-data[1];
}

*/