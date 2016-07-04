#include <iostream>
#include <math.h>
#include <complex>

#define npoints 10000000
#define niterate 10
#define PI 3.14159265359


	double w = 1.0;
	double k = 1.0;

double tex(const double & E, const double & x) //action integrando p(E,x)
{
	//if((2.0*E - w*w*x*x + (2.0/3.0)*k*x*x*x) < 0.0) std::cout << "errore " << 2.0*E - w*w*x*x + (2.0/3.0)*k*x*x*x <<std::endl <<std::flush;
	double temp = (2.0*E - w*w*x*x + (2.0/3.0)*k*x*x*x);
	return sqrt( temp<1e-10? 0.0 : 1.0/temp );
}

double dI_dE_simpson(double E, double a, double b, int N)
{
	// usare la funzione pex può restituire dei nan perché la funzione find root
	// approssima i valori numerici degli zeri di pex che è una sqrt quindi può
	// essere fuori dall'intervallo dove si hanno valori negativi per questo in
	// pex si controlla il valore prima di returnarlo
	double sum = 0.0;
	double dx= (b-a)/N;
		
	for (int i=0;i<N;i++)
	{
		sum += tex(E, a + i*dx) + 4.0*tex(E, a + 0.5*((i+1)*dx + i*dx) ) + tex(E, a + (i+1)*dx);
	}
		
	sum *= dx/6.0;
			
	return sum;
}

void find_roots(double & root1, double & root2, double & root3, const double a, const double b, const double c) //finds the first 2 roots of a-bx^2+cx^3
{
	std::complex<double> croot1, croot2, croot3, temp, ctemp, valueplus, valueminus;
		
	valueplus =  std::complex<double>(1.0,  sqrt(3.0));
	valueminus = std::complex<double>(1.0, -sqrt(3.0));
		
	ctemp =  std::sqrt( std::complex<double>( 27.0*a*a*pow(c, 16.0) - 4.0*a*b*b*b*pow(c, 14.0), 0.0 ));
	ctemp = 3.0*sqrt(3.0)*ctemp - 27.0*a*pow(c, 8.0) + 2.0*b*b*b*pow(c, 6.0);
	
	temp = std::pow(ctemp, 1.0/3.0);
	
	croot1 = b/(3.0*c) - valueminus*temp/(6.0*cbrt(2.0)*c*c*c) -  valueplus*b*b*c/(3.0*cbrt(4.0)*temp);
	croot2 = b/(3.0*c) -  valueplus*temp/(6.0*cbrt(2.0)*c*c*c) - valueminus*b*b*c/(3.0*cbrt(4.0)*temp);
	croot3 = b/(3.0*c) + ( temp/(cbrt(2.0)*c*c*c) + (cbrt(2.0)*b*b*c)/temp )/3.0;
	
	root1 = croot1.real();
	root2 = croot2.real();
	root3 = croot3.real();
}

double dI_dE_AGM(double a, double b, double c, int N) //ascending
{
	double alpha = asin(sqrt((b-a)/(c-a)));
	double x, y, z;
	x = 1.0;
	y = cos(alpha);
	z = sin(alpha);
		
	for(int i=0; i<N; i++) 
	{
		double x_temp = x;
		double y_temp = y;
		x = 0.5*(x_temp + y_temp);
		y = sqrt(x_temp * y_temp);
		z = 0.5*(x_temp - y_temp);
		
		if( z < 1e-10 ) break;
	}
	
	return PI/(2.0*x);
	//return phi/(x);
}

double dI_dE_landen(double a, double b, double c, int N) //ascending
{
	double alpha = asin(sqrt((b-a)/(c-a)));
	double k = 1.0;// + sin(alpha);
	
	for(int i=0; i<N; i++)
	{
		//phi = atan( cos(alpha)*tan(phi)) + phi;
		alpha = asin(2.0/(1.0 + cos(alpha)) - 1.0);
		
		k *= 1.0 + sin(alpha); 
		
		if( alpha < 1e-10 ) break;
		
	}
	return (PI*k/2.0)*2.0/sqrt(c-a)*sqrt(3.0/2.0);//pow(2.0, j);
}

double dI_dE_test(double a, double b, double c, int N) //ascending
{
	double alpha, phi, k, j;
	
	alpha = asin(sqrt((b-a)/(c-a)));
	phi = 3.0*PI/2.0;
	k = 1.0;
	
	for(int i=0; i<N; i++)
	{
		phi = atan( cos(alpha)*tan(phi)) + phi;
		alpha = asin(2.0/(1.0 + cos(alpha)) - 1.0);
		
		k *= 1.0 + sin(alpha); 
		
		if( alpha < 1e-5 ) {
			j = i;
			std::cout << "i " << i << std::endl;
			break;
		}
	}
	
	return (k*phi*2.0/sqrt(c-a))/pow(2.0, j);
}

int main(void)
{
	double Ei = 0.1;

	double root1, root2, root3;
		
	find_roots(root1, root2, root3, 2.0*Ei, w*w, 2.0*k/3.0);
	
	std::cout << root1 << " " << root2 << " " << root3 << std::endl;
	
	//double simpson = dI_dE_simpson(Ei, root1, root2, npoints);
	double landen  = dI_dE_landen(root1, root2, root3, niterate);
	double AGM = 2.0/sqrt(root3-root1)*sqrt(3.0/(2.0*k))*dI_dE_AGM(root1, root2, root3, niterate);

	
	//std::cout << "simpson = " << simpson << std::endl;
	std::cout << "landen  = " << landen  << std::endl;
	std::cout << "AGM  = " << AGM  << std::endl;
	return 0.0;
}