#include <vector>
#include <math.h>	// sqrt
#include <complex> //std::sqrt

#define PI 3.14159265359

double Klanden(double a, int N) //descending trasformation
{
	std::complex<double> alpha, k, temp;

	temp.real(a);
	temp.imag(0.0);

	alpha = std::asin(temp);
	k = 1.0;

	for(int i=0; i<N; i++)
	{
		alpha = std::asin(2.0/(1.0 + std::cos(alpha)) - 1.0);

		k *= 1.0 + std::sin(alpha); 

		if( std::abs(alpha) < 1e-10 ) break;
	}

	return PI*k.real()/2.0;
}

double Elanden(double a, int N) //descending trasformation
{
	std::complex<double> alpha, k, temp;

	temp.real(a);
	temp.imag(0.0);

	alpha = std::asin(temp);
	k = 1.0;
	std::vector<std::complex<double>> Nalpha;
	Nalpha.push_back(alpha);

	for(int i=0; i<N; i++)
	{
		alpha = std::asin(2.0/(1.0 + std::cos(alpha)) - 1.0);
		Nalpha.push_back(alpha);

		k *= 1.0 + std::sin(alpha); 

		if( std::abs(alpha) < 1e-10 ) break;		
	}

	std::complex<double>  c = 1.0 + 0.5*std::sin(Nalpha.back());
	for(int i=Nalpha.size()-2; i>0; i--) c = 1.0 + 0.5*std::sin(Nalpha[i])*c;

	std::complex<double> t = std::sin(Nalpha.front());
	c = 1.0 - 0.5*t*t*c;

	std::complex<double> res = PI*k*c/2.0;
	return res.real();
}