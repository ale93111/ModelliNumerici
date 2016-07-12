#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm> // std::max
#include <math.h>	// sqrt
#include <complex> //std::sqrt

#define npoints 1000
#define PI 3.14159265359

//#ifndef ENSEMBLE_HENON_H
//#define ENSEMBLE_HENON_H

#include"fourier.h"
#include"elliptic_int.h"
//#include"newton.h" included after definition of struct newton_func

std::random_device rd;
std::mt19937 generator(rd());
std::normal_distribution<double> gauss_distr(0.0, 1.0);
std::uniform_real_distribution<double> uniform_distr(-1, 1);

struct Newton_func
{
	//friend struct Ensemble;
	double I0, b, c;
	void (Newton_func::*fpointer)(double, double *, double *);
	
	void newtonfunc(double x, double *f, double *df);
	
	Newton_func(){}
	Newton_func(double I0i, double bi, double ci) : I0(I0i), b(bi), c(ci) {}
};
//include newton after definition of struct newton_func
#include"newton.h"

struct Ensemble
{
	//friend struct Newton_func;
		
	double *p, *q, *E, *action, *dI_dE;
	double w, k, epsilon, t;
	int Nparticles;
	
	double avg_energy()
	{
		double avg = 0.0;
		
		for(int i=0; i<Nparticles; i++) avg += E[i];
		
		return avg/Nparticles;
	}
	
	double avg_action()
	{
		double avg = 0.0;
		
		for(int i=0; i<Nparticles; i++) avg += action[i];
		
		return avg/Nparticles;
	}
	
	double avg_action_for_diffusion( const Ensemble & ensemble)
	{
		double avg = 0.0;
		
		for(int i=0; i<Nparticles; i++) avg += pow(action[i] - ensemble.action[i], 2.0);
		
		return avg/Nparticles;
	}
	
	double energy_max()
	{
		return w*w*w*w*w*w/(6.0*k*k);
	}
	
	void find_roots(double & root1, double & root2, double & root3, const double a, const double b, const double c) //finds the first 2 roots of a-bx^2+cx^3
	{
		std::complex<double> croot1, croot2, croot3, temp, ctemp, valueplus, valueminus;
		
		if( c == 0.0 ) 
		{
			root1 = -sqrt(a/b);
			root2 =  sqrt(a/b);
		}
		else 
		{
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
	}
	
	double energy(const double & q, const double & p)
	{
		return 0.5*p*p + 0.5*w*w*q*q - k*q*q*q/3.0;
	}
	
	double pex(const double & E, const double & x) //action integrando p(E,x)
	{
		double temp = 2.0*E - w*w*x*x + (2.0/3.0)*k*x*x*x;
		return sqrt( temp<1e-10? 0.0 : temp );
	}
	
	double tex(const double & E, const double & x) //action integrando p(E,x)
	{
		double temp = (2.0*E - w*w*x*x + (2.0/3.0)*k*x*x*x);
		return sqrt( temp<1e-10? 0.0 : 1.0/temp );
	}
	
	double action_elliptic(double root1, double root2, double root3)
	{
		return (2.0/(15.0*PI))*sqrt(2.0*k/3.0)*(sqrt(root2 - root1)*(root2 - root3)*(root1 - 2.0*root2 + root3)
										*Klanden(sqrt((root3 - root1)/(root2 - root1)), 10) 
								   + 2.0*sqrt(root2 - root1)*(root1*root1 + root2*root2 + root3*root3 - root2*root3 - root1*(root2 + root3))
										*Elanden(sqrt((root3 - root1)/(root2 - root1)), 10));
	}
	
	double dI_dE_elliptic(double root1, double root2, double root3)
	{
		return (1.0/PI)*sqrt(3.0/(2.0*k))*(2.0/sqrt(root3 - root1))*Klanden(sqrt((root2 - root1)/(root3 - root1)), 10);
	}
	
	void symplectic_advance(double & q, double & p, const double dt, const double t)
	{
		p = p - dt*( w*w*q - k*q*q );
		q = q + dt*p;
	}
		
	void advance( const double dt )
	{
		
		t +=dt;
		
		for(int i=0; i<Nparticles; i++)
		{
			symplectic_advance(q[i], p[i], dt/2.0, t);
		
			p[i] += epsilon*q[i]*sqrt(dt)*gauss_distr(generator);
		
			symplectic_advance(q[i], p[i], dt/2.0, t);
			
			E[i] = energy( q[i], p[i]); 
		
			double root1, root2, root3;
			
			find_roots(root1, root2, root3, 2.0*E[i], w*w, 2.0*k/3.0);		
			
			action[i] = action_elliptic(root1, root2, root3);
			dI_dE[i] = dI_dE_elliptic(root1, root2, root3);
		}
	}
	
	Ensemble(){}	
	void allocator( int N )
	{
		p = new double[N];
		q = new double[N];
		E = new double[N];
		action = new double[N];
		dI_dE = new double[N];	
	}
	
	Ensemble(int N, double wi, double ki, double epsiloni, double Ei) //: Ensemble( N )
	{
		allocator(N);
		Nparticles = N;
		w = wi;
		k = ki;
		epsilon = epsiloni;
		t = 0.0;
		
		double root1, root2, root3;
		
		find_roots(root1, root2, root3, 2.0*Ei, w*w, 2.0*k/3.0);
		
		//std::cout << "root 1 = " << root1 << " root2 = " << root2 << " root3 = " << root3 << std::endl;
		
		double actioni = action_elliptic(root1, root2, root3);
		double dI_dEi = dI_dE_elliptic(root1, root2, root3);
		
		for(int i=0; i<N; i++) 
		{
			E[i] = Ei;
			double nr = uniform_distr(generator);
			q[i] = (nr<0.0)? -nr*root1 : nr*root2;
			p[i] = pex(Ei, q[i]);
			action[i] = actioni;
			dI_dE[i] = dI_dEi;
		}
	}	
	
	Ensemble(int N, double wi, double ki, double epsiloni, double actioni, char Vi) //: Ensemble( N )
	{
		allocator(N);
		Nparticles = N;
		w = wi;
		k = ki;
		epsilon = epsiloni;
		t = 0.0;
		
		double root1, root2, root3;

		Newton_func fnew(actioni, w*w, 2.0*k/3.0);
	
		fnew.fpointer = &Newton_func::newtonfunc;
			
		double Ei = rtsafe(fnew, fnew.fpointer, 0.1*energy_max(), 0.9*energy_max(), 1e-5);
			
		find_roots(root1, root2, root3, 2.0*Ei, w*w, 2.0*k/3.0);
		
		double dI_dEi = dI_dE_elliptic(root1, root2, root3);

		
		for(int i=0; i<N; i++) 
		{
			E[i] = Ei;
			double nr = uniform_distr(generator);
			q[i] = (nr<0.0)? -nr*root1 : nr*root2;
			p[i] = pex(Ei, q[i]);
			action[i] = actioni;
			dI_dE[i] = dI_dEi;
		}
	}
	
	Ensemble& operator=(Ensemble&& other)
	{
		w = other.w;
		k = other.k;
		epsilon = other.epsilon;
		t = other.t;
		Nparticles = other.Nparticles;
		
		allocator(other.Nparticles);
		
		for(int i=0; i<Nparticles; i++) 
		{
			E[i] = other.E[i];
			q[i] = other.q[i];
			p[i] = other.p[i];
			action[i] = other.action[i];
			dI_dE[i] = other.dI_dE[i];
		}
		
		return *this;
	}
	
	Ensemble( const Ensemble& other ) //: Ensemble( other.Nparticles )
	{
		w = other.w;
		k = other.k;
		epsilon = other.epsilon;
		t = other.t;
		Nparticles = other.Nparticles;
		
		allocator(other.Nparticles);
		
		for(int i=0; i<Nparticles; i++) 
		{
			E[i] = other.E[i];
			q[i] = other.q[i];
			p[i] = other.p[i];
			action[i] = other.action[i];
			dI_dE[i] = other.dI_dE[i];
		}
	}
	
	~Ensemble()	
	{
		delete p, q, E, action, dI_dE;	
	}
};

void Newton_func::newtonfunc(double x, double *f, double *df)
{
	Ensemble ensemble;
	ensemble.allocator(1);
	ensemble.w = sqrt(b);
	ensemble.k = 3.0*c/2.0;
	
	double root1, root2, root3;
	ensemble.find_roots(root1, root2, root3, 2.0*x, b, c);
		
	f[0] = ensemble.action_elliptic(root1, root2, root3) - I0;
	df[0] = ensemble.dI_dE_elliptic(root1, root2, root3);
}

//#endif