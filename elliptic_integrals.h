#include <vector>
#include <math.h>	// sqrt
#include <complex> //std::sqrt 

#define PI 3.14159265359

//descending landen trasformations, see abramowitz and stegun
//generalized to complex values of the modular angle, they return only the real part of the result

//IMPORTANT: the argument of these functions is the modular angle
//check wolframalpha page on elliptic integrals or abramotiz and stegun for possible other conventions

//complete elliptic integral of the first kind
double Klanden(double a, int N) 
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

//complete elliptic integral of the second kind
double Elanden(double a, int N) 
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

