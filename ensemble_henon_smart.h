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
//#include"newton.h" 

std::random_device rd;
std::mt19937 generator(rd());
std::normal_distribution<double> noise_gauss_distr(0.0, 1.0);
std::uniform_real_distribution<double> uniform_distr(-1, 1);

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

//void newtonfunc(double x, double *f, double *df, double I0);

struct Ensemble
{
	//friend struct Newton_func;
		
	//double *p, *q, *E, *action, *dI_dE;
	std::vector<double> p, q, E, action, dI_dE;
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
	
	double avg_dI_dE()
	{
		double avg = 0.0;
		
		for(int i=0; i<Nparticles; i++) avg += dI_dE[i];
		
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
	
	void check_energy(int index)
	{
		// è poco efficiente eliminare valori che non sono allocati alla fine!
		if( E[index] < 0.1*energy_max() || E[index] > 0.9*energy_max() )
		{
			p.erase( p.begin() + index );
			q.erase( q.begin() + index );
			E.erase( E.begin() + index );
			action.erase( action.begin() + index );
			dI_dE.erase( dI_dE.begin() + index );
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
		
			p[i] += epsilon*q[i]*sqrt(dt)*noise_gauss_distr(generator);
		
			symplectic_advance(q[i], p[i], dt/2.0, t);
			
			E[i] = energy( q[i], p[i]); 
		
			double root1, root2, root3;
			
			find_roots(root1, root2, root3, 2.0*E[i], w*w, 2.0*k/3.0);		
			
			action[i] = action_elliptic(root1, root2, root3);
			dI_dE[i] = dI_dE_elliptic(root1, root2, root3);
		}
	}
	
	void newtonfunc(double x, double *f, double *df, double I0)
	{	
		double root1, root2, root3;
		find_roots(root1, root2, root3, 2.0*x, w*w, 2.0*k/3.0);
		
		f[0] = action_elliptic(root1, root2, root3) - I0;
		df[0] = dI_dE_elliptic(root1, root2, root3);
	}
	
	double rtsafe(double x1, double x2, double xacc, double offset);
	
	void allocator( int N )
	{
		p.resize( N, 0);
		q.resize( N, 0);
		E.resize( N, 0);
		action.resize( N, 0);
		dI_dE.resize( N, 0);
	}
	
	Ensemble(){}
	
	Ensemble(int N, double wi, double ki, double epsiloni, double Ei) : Nparticles(N), w(wi), k(ki), epsilon(epsiloni) 
	{
		allocator(N);
		t = 0.0;
		
		double root1, root2, root3;
		
		find_roots(root1, root2, root3, 2.0*Ei, w*w, 2.0*k/3.0);
				
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
	
	Ensemble(int N, double wi, double ki, double epsiloni, double actioni, char Vi) : Nparticles(N), w(wi), k(ki), epsilon(epsiloni)
	{
		allocator(N);
		t = 0.0;
		
		double root1, root2, root3;
			
		double Ei = rtsafe( 0.1*energy_max(), 0.9*energy_max(), 1e-5, actioni);
			
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
	
	Ensemble(int N, double wi, double ki, double epsiloni, double mean, double devstd, char Vi) : Nparticles(N), w(wi), k(ki), epsilon(epsiloni)
	{
		allocator(N);
		t = 0.0;
		
		std::normal_distribution<double> gauss_distr(mean, devstd);
		
		double Ea = 0.1*energy_max();
		double Eb = 0.9*energy_max();
		
		for(int i=0; i<N; i++)
		{
			double Ei, actioni;
			double root1, root2, root3;
			
			do
			{
				actioni = gauss_distr(generator);
		
				Ei = rtsafe(Ea, Eb, 1e-5, actioni);
			} 
			while( Ei < Ea || Ei > Eb);
			
			find_roots(root1, root2, root3, 2.0*Ei, w*w, 2.0*k/3.0);
			
			dI_dE[i] = dI_dE_elliptic(root1, root2, root3);
			
			E[i] = Ei;
			double nr = uniform_distr(generator);
			q[i] = (nr<0.0)? -nr*root1 : nr*root2;
			p[i] = pex(Ei, q[i]);
			action[i] = actioni;
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
	
	Ensemble( const Ensemble& other )
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
	
	//~Ensemble()	{ delete p, q, E, action, dI_dE;	}
};

//from Numerical recipe in c
double Ensemble::rtsafe(double x1, double x2, double xacc, double offset = 0.0)
{
	//Maximum allowed number of iterations.
	const int MAXIT = 100;
	int j;
	double df,dx,dxold,f,fh,fl;
	double temp,xh,xl,rts;
	
	newtonfunc(x1,&fl,&df,offset);
	newtonfunc(x2,&fh,&df,offset);
	
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
		std::cout << "Root must be bracketed in rtsafe" << std::endl;
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		//Orient the search so that f(xl)<0. 
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	}
	rts=0.5*(x1+x2);
	//Initialize the  guess  for root,
	dxold=fabs(x2-x1);
	//the “stepsize  before  last,”
	dx=dxold;
	//and  the  last  step.
	newtonfunc(rts,&f,&df,offset);
	for (j=1;j<=MAXIT;j++) {
		//Loop  over  allowed  iterations.
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
		//Bisect if Newton out of range,
		|| (fabs(2.0*f) > fabs(dxold*df))) {
		//or not decreasing fast enough.
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			if (xl == rts) return rts;
			//Change  in root  is negligible.
		} else {
			//Newton  step acceptable.  Take it.
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts;
		//Convergence  criterion.
		newtonfunc(rts,&f,&df,offset);
		//The  one  new  function  evaluation  per  iteration.
		if (f < 0.0)
		//Maintain the  bracket  on the root.
			xl=rts;
		else
			xh=rts;
	}
	std::cout << "Maximum number of iterations exceeded in rtsafe" << std::endl;
	return 0.0;
	//Never  get  here.
}

//#endif