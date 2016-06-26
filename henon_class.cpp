#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm> // std::max
#include <math.h>	// sqrt
#include <complex> //std::sqrt

#include"fourier.h"
#include"newton.h"

#define npoints 1000
#define PI 3.14159265359

std::random_device rd;
std::mt19937 generator(rd());
std::normal_distribution<double> gauss_distr(0.0, 1.0);
std::uniform_real_distribution<double> uniform_distr(-1, 1);


void pol(double x, double *func, double *dfunc, const double a, const double b, const double c)
{
	*func  = a - b*x*x + c*x*x*x;
	*dfunc =  - 2.0*b*x + 3.0*c*x*x;
}


struct Ensemble
{
	double *p, *q, *E, *action, *dI_dE;
	double w, k, epsilon, t;
	int Nparticles;
	void (*pfpolinomio)(double, double *, double *, const double , const double , const double ); 	
	
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
		
		//root1 = rtsafe(pfpolinomio, root1 - 0.01*root1, root1 + 0.01*root1, 1e-6, a, b, c);
		//root2 = rtsafe(pfpolinomio, root2 - 0.01*root2, root2 + 0.01*root2, 1e-6, a, b, c);
		
		std::cout << "root1 = " << root1 << "  root2 = " << root2 << std::endl;
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
		double temp = (2.0*E - w*w*x*x + (2.0/3.0)*k*x*x*x);
		return sqrt( temp<1e-10? 0.0 : 1.0/temp );
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
			
		return sum/PI;
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
			
		return sum/PI;
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
			dI_dE[i] = dI_dE_simpson(E[i], root1, root2, npoints);
		}
	}
	
	Ensemble(){}	
	Ensemble( int N )
	{
		p = new double[N];
		q = new double[N];
		E = new double[N];
		action = new double[N];
		dI_dE = new double[N];	
	}
	
	Ensemble(int N, double wi, double ki, double epsiloni, double Ei)
	{
		pfpolinomio = pol;
		Nparticles = N;
		w = wi;
		k = ki;
		epsilon = epsiloni;
		
		p = new double[N];
		q = new double[N];
		E = new double[N];
		action = new double[N];
		dI_dE = new double[N];
		
		t = 0.0;
		
		double root1, root2;
		
		find_roots(root1, root2, 2.0*Ei, w*w, 2.0*k/3.0);
				std::cout << "ciao sono il costruttore" << std::endl;
		double actioni = action_simpson(Ei, root1, root2, npoints);
		double dI_dEi = dI_dE_simpson(Ei, root1, root2, npoints);
		
		for(int i=0; i<N; i++) 
		{
			E[i] = Ei;
			q[i] = 0.0;
			p[i] = pex(Ei, q[i]);
			action[i] = actioni;
			dI_dE[i] = dI_dEi;
		}
	}	
	
	Ensemble( const Ensemble& other ) : Ensemble( other.Nparticles )
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
			dI_dE[i] = other.dI_dE[i];
		}
	}
	
	~Ensemble()	
	{
		delete p, q, E, action, dI_dE;	
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

	int Nensemble = 1;
	int nsteps = 10;
	int Ndynamic = 512;
	
	double w = 1.0;
	double k = 1.0;
	double dt = 0.01;
	double dtdynamic = 0.001;
	double epsilon = 0.1;
	
	double E = 0.1;
	double Emax = w*w*w*w*w*w/(6.0*k*k);
	
	std::cout << "Emax = " << Emax << std::endl;
	std::cout << std::endl;

	double coeff_drift, coeff_diffusion, coeff_diffusion_theoretical;
	std::vector<double> Ncoeff_diffusion, Ntime;
	std::vector<double> Ncoeff_diffusion_dynamic;
	
	double *q_p_f;
	q_p_f = new double[Ndynamic*Nensemble];
	
	
	Ensemble ensemble(Nensemble, w, k, epsilon, E);
	
	Ensemble ensemble_0 = ensemble;
	
	std::cout << "t iniziale = " << ensemble.t << "\t" << "E[0] iniziale = " << ensemble.E[0] << std::endl
			  << "azione[0] iniziale = " << ensemble.action[0] << "\t" << "dI_dE[0] iniziale = " << ensemble.dI_dE[0] << std::endl;
	
	double avg_action_0 = ensemble.avg_action();
	
	for(int i=0; i<nsteps; i++)
	{
		ensemble.advance(dt);
		
		Ntime.push_back(ensemble.t);
		Ncoeff_diffusion.push_back( ensemble.avg_action_for_diffusion(ensemble_0)/dt ); 
		
		Ensemble ensemble_temp = ensemble;
		//propaga la dinamica simplettica per Ndynamic passi
		for(int j=0; j<Ndynamic; j++)
		{
			//per ogni particella
			for(int n=0; n<ensemble_temp.Nparticles; n++)
			{
				ensemble_temp.symplectic_advance( ensemble_temp.q[n], ensemble_temp.p[n], dtdynamic, ensemble_temp.t);
				
				//costruisco un segnale di cui vado a fare la FFT
				q_p_f[n*Ndynamic + j] = pow( ensemble_temp.q[n]*ensemble_temp.p[n]/ensemble_temp.dI_dE[n] , 2.0); // ;
			}
			
			ensemble_temp.t += dtdynamic;
		}
		
		double avg_diffusion_temp = 0.0;
		
		for(int n=0; n<ensemble_temp.Nparticles; n++)
		{
			//FFT per la traiettoria di ogni particella
			realft( &q_p_f[n*Ndynamic], Ndynamic );
			
			avg_diffusion_temp += q_p_f[n*Ndynamic];	
			//avg_diffusion_temp += sqrt(q_p_f[n*Ndynamic]*q_p_f[n*Ndynamic] + q_p_f[n*Ndynamic + 1]*q_p_f[n*Ndynamic + 1]) / Ndynamic;	
		}
		avg_diffusion_temp /= Ndynamic*ensemble_temp.Nparticles;
		Ncoeff_diffusion_dynamic.push_back( epsilon*epsilon*avg_diffusion_temp );		
	}
	
	coeff_drift = (ensemble.avg_action() - avg_action_0)/dt;
	coeff_diffusion = slope_linear_regression(Ntime, Ncoeff_diffusion);
	//coeff_diffusion_theoretical = slope_linear_regression(Ntime, Ncoeff_diffusion_dynamic);
	coeff_diffusion_theoretical = accumulate(Ncoeff_diffusion_dynamic.begin(), Ncoeff_diffusion_dynamic.end(), 0.0) / Ncoeff_diffusion_dynamic.size();
		
		
	std::cout << "t finale = " << ensemble.t << "\t" << "E[0] finale = " << ensemble.E[0] << std::endl
			  << "azione[0] finale = " << ensemble.action[0] << "\t" << "dI_dE[0] finale = " << ensemble.dI_dE[0]<< std::endl;
	
	
	
	std::cout << std::endl << std::endl;
	std::cout << "coeff_drift = " << coeff_drift << std::endl;
	std::cout << "coeff_diffusion = " << coeff_diffusion << std::endl;
	std::cout << "coeff_diffusion_theoretical = " << coeff_diffusion_theoretical << std::endl;
	
	//for(int i=0; i<Ncoeff_diffusion.size(); i++) 		 std::cout << Ncoeff_diffusion[i] 		  << std::endl;
	std::cout << std::endl;
	//for(int i=0; i<Ncoeff_diffusion_dynamic.size(); i++) std::cout << Ncoeff_diffusion_dynamic[i] << std::endl;
	
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