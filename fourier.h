#include <math.h>	// sqrt

#define PI 3.14159265359
#define TWOPI	(2.0*PI)

//from Numerical recipe in c, computes the fourier transform of data (of size nn, which should be a power of 2)
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

//zero padding of complex components to compute the fourier transform of a real signal
void realft(double* data, unsigned long n)
{
	double *X = new double[2*n+1]; 
	for(int i=0; i<n; i++)
	{
		X[2*i+1] = data[i];
		X[2*i+2] = 0.0;
	}	
	
	four1(X, n, 1);
	
	for(int i=0; i<n; i++) data[i] = X[i+1];
}

