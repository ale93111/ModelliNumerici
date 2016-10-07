#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm> // std::max
#include <math.h>	// sqrt
#include <complex> //std::sqrt
#include <omp.h> //OpenMP 

#include"fourier.h"
#include"elliptic_integrals.h"

#define PI 3.14159265359

std::random_device rd;
//std::mt19937 generator(rd());

#ifdef _OPENMP
std::mt19937 generator(omp_get_thread_num());
#else
std::mt19937 generator(0);//rd());
#endif

std::normal_distribution<double> noise_gauss_distr(0.0, 1.0);
std::uniform_real_distribution<double> uniform_distr(-1, 1);

//finds the real part of the roots of the polynomial a-bx^2+cx^3
void find_roots(double & root1, double & root2, double & root3, const double a, const double b, const double c) 
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


struct Ensemble
{
	std::vector<double> p, q, E, action, dI_dE;
	//store a deleted list of the indexes of the particles that have energies outside of the defined range of the problem (10% ~ 90% of the max energy)
	std::vector<int> deletedlist;
	double w, k, epsilon, t;
	int Nparticles;
	
	//compute the average energy of the ensemble of particles
	double avg_energy()
	{
		double avg = 0.0;
		
		//#pragma omp parallel for reduction(+:avg)
		for(int i=0; i<Nparticles; i++) avg += E[i];
		
		return avg/Nparticles;
	}
	
	//compute the average action of the ensemble of particles
	double avg_action()
	{
		double avg = 0.0;
		
		//#pragma omp parallel for reduction(+:avg)
		for(int i=0; i<Nparticles; i++) avg += action[i];
		
		return avg/Nparticles;
	}
	
	//compute the average dI/dE of the ensemble of particles
	double avg_dI_dE()
	{
		double avg = 0.0;
		
		//#pragma omp parallel for reduction(+:avg)
		for(int i=0; i<Nparticles; i++) avg += dI_dE[i];
		
		return avg/Nparticles;
	}
	
	//compute the second momentum of the action distribution between two ensembles (used to find the diffusion coefficient) 
	double avg_action_for_diffusion( const Ensemble & ensemble)
	{
		double avg = 0.0;
		
		//#pragma omp parallel for reduction(+:avg)
		for(int i=0; i<Nparticles; i++) avg += pow(action[i] - ensemble.action[i], 2.0);
		
		return avg/Nparticles;
	}
	
	//compute the max energy that corresponds to the stability region of the phase space (boundary=separatrix)
	double energy_max()
	{
		return w*w*w*w*w*w/(6.0*k*k);
	}
	
	//check the value of the energy of the particle with index=index, if it's outside the range of allowed energies deletes the particle from the ensemble
	void check_energy(int index)
	{
		//check the range of allowed energies
		if( E[index] < 0.1*energy_max() || E[index] > 0.9*energy_max() )
		{
			if(p.empty()) 
			{
				Nparticles = 0;
				return;
			}
			
			//erase particle
			p.erase( p.begin() + index );
			q.erase( q.begin() + index );
			E.erase( E.begin() + index );
			action.erase( action.begin() + index );
			dI_dE.erase( dI_dE.begin() + index );
		
			//update number of particles
			Nparticles = p.size();
			//store index
			deletedlist.push_back(index);
		}
	}

	//compute the energy of a particle from p and q
	double energy(const double & q, const double & p)
	{
		return 0.5*p*p + 0.5*w*w*q*q - k*q*q*q/3.0;
	}
	
	//action integrand p(E,x)
	double pex(const double & E, const double & x)
	{
		double temp = 2.0*E - w*w*x*x + (2.0/3.0)*k*x*x*x;
		return sqrt( temp<1e-10? 0.0 : temp ); //fix to avoid numerical troubles
	}
	
	double tex(const double & E, const double & x) // t(E,x)
	{
		double temp = (2.0*E - w*w*x*x + (2.0/3.0)*k*x*x*x);
		return sqrt( temp<1e-10? 0.0 : 1.0/temp );
	}
	
	//compute the analytical action integral using complete elliptic integrals
	double action_elliptic(double root1, double root2, double root3)
	{
		return (2.0/(15.0*PI))*sqrt(2.0*k/3.0)*(sqrt(root2 - root1)*(root2 - root3)*(root1 - 2.0*root2 + root3)
										*Klanden(sqrt((root3 - root1)/(root2 - root1)), 10) 
			  			 + 2.0*sqrt(root2 - root1)*(root1*root1 + root2*root2 + root3*root3 - root2*root3 - root1*(root2 + root3))
										*Elanden(sqrt((root3 - root1)/(root2 - root1)), 10));
	}
	
	//compute the analytical dI/dE integral using complete elliptic integrals
	double dI_dE_elliptic(double root1, double root2, double root3)
	{
		return (1.0/PI)*sqrt(3.0/(2.0*k))*(2.0/sqrt(root3 - root1))*Klanden(sqrt((root2 - root1)/(root3 - root1)), 10);
	}
	
	//symplectic integrator of the motion equations (splitting method)
	void symplectic_advance(double & q, double & p, const double dt, const double t)
	{
		p = p - dt*( w*w*q - k*q*q );
		q = q + dt*p;
	}
	
	//algorithm to advance the dynamics with noise (see report)
	void advance( const double dt )
	{
		
		t +=dt;
		
		#pragma omp parallel for 
		for(int i=0; i<Nparticles; i++)
		{
			symplectic_advance(q[i], p[i], dt/2.0, t);
			
			p[i] += epsilon*q[i]*sqrt(dt)*noise_gauss_distr(generator);
		
			symplectic_advance(q[i], p[i], dt/2.0, t);
			
			//find the new energy, action and dI/dE
			E[i] = energy( q[i], p[i]); 
		
			double root1, root2, root3;
			
			find_roots(root1, root2, root3, 2.0*E[i], w*w, 2.0*k/3.0);		
			
			action[i] = action_elliptic(root1, root2, root3);
			dI_dE[i] = dI_dE_elliptic(root1, root2, root3);
			
			//check energy of particles and delete the particles outside range
			check_energy(i);
		}
	}
	
	//this is the function used by the newton method, which is used to initialize the ensemble with an action value or distribution.
	//We know how to compute analitically the action I of a particle from its energy E, but not viceversa.
	//So we know I(E) but not E(I). To find E(I) we use the newton method to find the zero of the function I(E)-I0
	//where I0 is the desired action value that we want to use to initialize the values of the particles
	void newtonfunc(double x, double *f, double *df, double I0)
	{	
		double root1, root2, root3;
		find_roots(root1, root2, root3, 2.0*x, w*w, 2.0*k/3.0);
		
		f[0] = action_elliptic(root1, root2, root3) - I0;
		df[0] = dI_dE_elliptic(root1, root2, root3);
	}
	
	//from Numerical recipe in c, the newton method applied to newtonfunc defined above
	double rtsafe(double x1, double x2, double xacc, double offset);
	//it's defined below
	
	void allocator( int N )
	{
		p.resize( N, 0);
		q.resize( N, 0);
		E.resize( N, 0);
		action.resize( N, 0);
		dI_dE.resize( N, 0);
	}
	
	Ensemble(){}
	
	//constructor used to initialize the ensemble with the same value of energy for all the particles
	Ensemble(int N, double wi, double ki, double epsiloni, double Ei) : Nparticles(N), w(wi), k(ki), epsilon(epsiloni) 
	{
		allocator(N);
		t = 0.0;
		
		double root1, root2, root3;
		
		find_roots(root1, root2, root3, 2.0*Ei, w*w, 2.0*k/3.0);
				
		double actioni = action_elliptic(root1, root2, root3);
		double dI_dEi = dI_dE_elliptic(root1, root2, root3);
		
		#pragma omp parallel for 
		for(int i=0; i<N; i++) 
		{
			E[i] = Ei;
			
			double nr = uniform_distr(generator);
			
			q[i] = (nr<0.0)? -nr*root1 : nr*root2;
			
			nr = uniform_distr(generator);
			
			p[i] = ( (nr<0.0)? -1.0 : 1.0 )*pex(Ei, q[i]);
			action[i] = actioni;
			dI_dE[i] = dI_dEi;
		}
	}	
	
	//constructor overload to initialize the ensemble with the same value of action for all the particles
	//overload: add a char as last argument
	Ensemble(int N, double wi, double ki, double epsiloni, double actioni, char Vi) : Nparticles(N), w(wi), k(ki), epsilon(epsiloni)
	{
		allocator(N);
		t = 0.0;
		
		double root1, root2, root3;
			
		double Ei = rtsafe( 0.1*energy_max(), 0.9*energy_max(), 1e-5, actioni);
			
		find_roots(root1, root2, root3, 2.0*Ei, w*w, 2.0*k/3.0);
		
		double dI_dEi = dI_dE_elliptic(root1, root2, root3);

		#pragma omp parallel for 
		for(int i=0; i<N; i++) 
		{
			E[i] = Ei;
			
			double nr = uniform_distr(generator);
			
			q[i] = (nr<0.0)? -nr*root1 : nr*root2;
			
			nr = uniform_distr(generator);
			
			p[i] = ( (nr<0.0)? -1.0 : 1.0 )*pex(Ei, q[i]);
			
			action[i] = actioni;
			dI_dE[i] = dI_dEi;
		}
	}
	
	//overload constructor to initialize the ensemble with a Gaussian action distribution  
	Ensemble(int N, double wi, double ki, double epsiloni, double mean, double devstd, char Vi) : Nparticles(N), w(wi), k(ki), epsilon(epsiloni)
	{
		allocator(N);
		t = 0.0;
		
		std::normal_distribution<double> gauss_distr(mean, devstd);
		
		double Ea = 0.1*energy_max();
		double Eb = 0.9*energy_max();
		
		#pragma omp parallel for 
		for(int i=0; i<N; i++)
		{
			double Ei, actioni;
			double root1, root2, root3;
			
			//exclude values of action sampled from the distribution that correspond to values of energy outside range
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
			
			nr = uniform_distr(generator);
			
			p[i] = ( (nr<0.0)? -1.0 : 1.0 )*pex(Ei, q[i]);
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
	
	Ensemble& operator=(Ensemble& other)
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
};

//from Numerical recipe in c, the newton methon applied to newtonfunc
//[x1,x2] is the interval where the algorithm looks for the zero of the function
//xacc is the desired accuracy of the solution
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
