#include <math.h>	// sqrt

void four1(double* data, unsigned long nn)
{
	unsigned long n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;
 
	// reverse-binary reindexing
	n = nn<<1;
	j=1;
	for (i=1; i<n; i+=2) {
		if (j>i) {
			std::swap(data[j-1], data[i-1]);
			std::swap(data[j], data[i]);
		}
		m = nn;
		while (m>=2 && j>m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	};
 
	// here begins the Danielson-Lanczos section
	mmax=2;
	while (n>mmax) {
		istep = mmax<<1;
		theta = -(2*M_PI/mmax);
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m=1; m < mmax; m += 2) {
			for (i=m; i <= n; i += istep) {
				j=i+mmax;
				tempr = wr*data[j-1] - wi*data[j];
				tempi = wr * data[j] + wi*data[j-1];
 
				data[j-1] = data[i-1] - tempr;
				data[j] = data[i] - tempi;
				data[i-1] += tempr;
				data[i] += tempi;
			}
			wtemp=wr;
			wr += wr*wpr - wi*wpi;
			wi += wi*wpr + wtemp*wpi;
		}
		mmax=istep;
	}
}

//Calculates the Fourier transform of a set of n real-valued data points. 
//Replaces this data (which is stored in array data[1..n]) 
//by the positive frequency half of its complex Fourier transform.
//The real-valued first and last components of the complex transform are returned as elements
//data[1] and data[2], respectively. n must be a power of 2. This routine also calculates the
//inverse transform of a complex data array if it is the transform of real data. 
//(Result in this case must be multiplied by 2/n.)
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