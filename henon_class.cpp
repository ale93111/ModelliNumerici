#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm> // std::max
#include <math.h>    // sqrt
#include <complex> //std::sqrt

#define npoints 1000
#define PI 3.14159265359

std::random_device rd;
std::mt19937 generator(rd());
std::normal_distribution<double> gauss_distr(0.0, 1.0);
std::uniform_real_distribution<double> uniform_distr(-1, 1);


struct Ensemble
{
	double *p, *q, *E, *action, *period;
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
	
	void find_roots(double & root1, double & root2, const double a, const double b, const double c) //finds the first 2 roots of a-bx^2+cx^3
	{
		std::complex<double> croot1, croot2, temp, ctemp, valueplus, valueminus;
		
		valueplus =  std::complex<double>(1.0,  sqrt(3.0));
		valueminus = std::complex<double>(1.0, -sqrt(3.0));
		
		ctemp =  std::sqrt( std::complex<double>( 27.0*a*a*pow(c, 16.0) - 4.0*a*b*b*b*pow(c, 14.0), 0.0 ));
		ctemp = 3.0*sqrt(3.0)*ctemp - 27.0*a*pow(c, 8.0) + 2.0*b*b*b*pow(c, 6.0);
	
		temp = std::pow(ctemp, 1.0/3.0);
	
		croot1 = b/(3.0*c) - valueminus*temp/(6.0*cbrt(2.0)*c*c*c) -  valueplus*b*b*c/(3.0*cbrt(4.0)*temp);
		croot2 = b/(3.0*c) -  valueplus*temp/(6.0*cbrt(2.0)*c*c*c) - valueminus*b*b*c/(3.0*cbrt(4.0)*temp);
	
		root1 = croot1.real();
		root2 = croot2.real();
	}
	
	
	double energy(const double & q, const double & p)
	{
		return 0.5*p*p + 0.5*w*w*q*q - k*q*q*q/3.0;
	}
	
	double pex(const double & E, const double & x) //action integrando p(E,x)
	{
		//if((2.0*E - w*w*x*x + (2.0/3.0)*k*x*x*x) < 0.0) std::cout << "errore " << 2.0*E - w*w*x*x + (2.0/3.0)*k*x*x*x <<std::endl <<std::flush;
		double temp = 2.0*E - w*w*x*x + (2.0/3.0)*k*x*x*x;
		return sqrt( temp<1e-10? 0.0 : temp );
	}
	
	double tex(const double & E, const double & x) //action integrando p(E,x)
	{
		//if((2.0*E - w*w*x*x + (2.0/3.0)*k*x*x*x) < 0.0) std::cout << "errore " << 2.0*E - w*w*x*x + (2.0/3.0)*k*x*x*x <<std::endl <<std::flush;
		double temp = 1.0/(2.0*E - w*w*x*x + (2.0/3.0)*k*x*x*x);
		return sqrt( temp<1e-10? 0.0 : temp );
	}
	
	double action_simpson(double E, double a, double b, int N)
	{
		// usare la funzione pex può restituire dei nan perché la funzione find root
		// approssima i valori numerici degli zeri di pex che è una sqrt quindi può
		// essere fuori dall'intervallo dove si hanno valori negativi per questo in
		// pex si controlla il valore prima di returnarlo
		double sum = 0.0;
		double dx= (b-a)/N;
		
		for (int i=0;i<N;i++)
		{
			sum += pex(E, a + i*dx) + 4.0*pex(E, a + 0.5*((i+1)*dx + i*dx) ) + pex(E, a + (i+1)*dx);
		}
		
		sum *= dx/6.0;
			
		return 2.0*sum;
	}
	
	double period_simpson(double E, double a, double b, int N)
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
			
		return 2.0*sum;
	}
	
	void symplectic_advance(double & q, double & p, const double dt, const double t)
	{
		p = p + dt*( w*w*q - k*q*q );
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
			
			E[i] = energy( q[i], p[i]); //NON è il contrario!!
		
			double root1, root2;
			
			find_roots(root1, root2, 2.0*E[i], w*w, 2.0*k/3.0);		
			
			action[i] = action_simpson(E[i], root1, root2, npoints);
			//period[i] = period_simpson(E[i], root1, root2, npoints);
		}
	}
	
		
	Ensemble(){}
	Ensemble(int N, double wi, double ki, double epsiloni, double Ei)
	{
		Nparticles = N;
		w = wi;
		k = ki;
		epsilon = epsiloni;
		
		p = new double[N];
		q = new double[N];
		E = new double[N];
		action = new double[N];
		period = new double[N];
		
		t = 0.0;
		
		double root1, root2;
		
		find_roots(root1, root2, 2.0*Ei, w*w, 2.0*k/3.0);

		double actioni = action_simpson(Ei, root1, root2, npoints);
		double periodi = period_simpson(Ei, root1, root2, npoints);
		
		for(int i=0; i<N; i++) 
		{
			E[i] = Ei;
			q[i] = 0.0;
			p[i] = pex(Ei, q[i]);
			action[i] = actioni;
			period[i] = periodi;
		}
	}	
	
	Ensemble( const Ensemble& other ) : Ensemble(other.Nparticles, other.w, other.k, other.epsilon, other.E[0])
	{
		w = other.w;
		k = other.k;
		epsilon = other.epsilon;
		t = other.t;
		Nparticles = other.Nparticles;
		
		for(int i=0; i<Nparticles; i++) 
		{
			E[i] = other.E[i];
			q[i] = other.q[i];
			p[i] = other.p[i];
			action[i] = other.action[i];
			period[i] = other.period[i];
		}
	}
	
	~Ensemble()	
	{
		delete p, q, E, action, period;	
	}
};


double slope_linear_regression(const std::vector<double>& x, const std::vector<double>& y)
{
	int n = x.size();
	
	double avgX = accumulate(x.begin(), x.end(), 0.0) / n;
    double avgY = accumulate(y.begin(), y.end(), 0.0) / n;
	
	double Sxy = 0.0;
    double Sxx = 0.0;

    for(int i=0; i<n; ++i){
        Sxy += (x[i] - avgX) * (y[i] - avgY);
        Sxx += (x[i] - avgX) * (x[i] - avgX);
    }

    return Sxy / Sxx;
}

int main(int argc, char *argv[])
{

	int Nensemble = 1000;
	int nsteps = 10;
	
	double w = 1.0;
	double k = 1.0;
	double dt = 0.01;
	double epsilon = 0.1;
	
	double E = 0.1;
	double Emax = w*w*w*w*w*w/(6.0*k*k);
	
	std::cout << "Emax = " << Emax << std::endl;
	std::cout << std::endl;

	double coeff_drift, coeff_diffusion;
	std::vector<double> Ncoeff_diffusion, Ntime;
	
	Ensemble ensemble(Nensemble, w, k, epsilon, E);
	
	Ensemble ensemble_0 = ensemble;
	
	std::cout << "t iniziale = " << ensemble.t << "\t" << "E[0] iniziale = " << ensemble.E[0] << std::endl
			  << "azione[0] iniziale = " << ensemble.action[0] << "\t" << "period[0] iniziale = " << ensemble.period[0] << std::endl;
	
	double avg_action_0 = ensemble.avg_action();
	
	for(int i=0; i<nsteps; i++)
	{
		ensemble.advance(dt);
		
		Ntime.push_back(ensemble.t);
		Ncoeff_diffusion.push_back( ensemble.avg_action_for_diffusion(ensemble_0)/dt ); 
		//Ncoeff_diffusion.push_back( pow( ensemble.avg_action() - avg_action_0 , 2.0)/dt ); //check formula of average values
	}
	
	coeff_drift = (ensemble.avg_action() - avg_action_0)/dt;
	coeff_diffusion = slope_linear_regression(Ntime, Ncoeff_diffusion);
	
	
	std::cout << "t finale = " << ensemble.t << "\t" << "E[0] finale = " << ensemble.E[0] << std::endl
			  << "azione[0] finale = " << ensemble.action[0] << "\t" << "period[0] finale = " << ensemble.period[0]<< std::endl;
	
	
	
	std::cout << std::endl << std::endl;
	std::cout << "coeff_drift = " << coeff_drift << std::endl;
	std::cout << "coeff_diffusion = " << coeff_diffusion << std::endl;
	
	//for(int i=0; i<Ncoeff_diffusion.size(); i++) std::cout << Ncoeff_diffusion[i] << std::endl;
	
	
	/*
	std::ofstream output;	
	output.open("out.txt");
	
	for(int i=0; i<Nensemble; i++)
	{
		output << ensemble_0.action[i] << "\t" << ensemble.action[i] << std::endl; 	
	}
	
	output.close();
	*/
	
	std::cout << std::endl;
	std::cout << "Fatto!" << std::endl;
	
	return 0;
	
}