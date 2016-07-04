#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm> // std::max
#include <math.h>	// sqrt
#include <complex> //std::sqrt

#include"fourier.h"
//#include"newton.h"

#define npoints 1000
#define niterate 10
#define PI 3.14159265359

std::random_device rd;
std::mt19937 generator(rd());
std::normal_distribution<double> gauss_distr(0.0, 1.0);
std::uniform_real_distribution<double> uniform_distr(-1, 1);

struct Ensemble
{
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
	
	double Klanden(double a, int N) //descending
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
	double Elanden(double a, int N) //descending
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
		
		double root1, root2, root3;
		
		find_roots(root1, root2, root3, 2.0*Ei, w*w, 2.0*k/3.0);
		
		std::cout << "root 1 = " << root1 << " " << "root2 = " << root2 << "root3 = " << root3 << std::endl;
		
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

double theoretical_diffusion(Ensemble ensemble_temp, int Ndynamic, double *q_p_f)
{	
	//propaga la dinamica simplettica per Ndynamic passi
	for(int j=0; j<Ndynamic; j++)
	{		
		//per ogni particella
		for(int n=0; n<ensemble_temp.Nparticles; n++)
		{
			ensemble_temp.symplectic_advance( ensemble_temp.q[n], ensemble_temp.p[n], 2.0*PI*ensemble_temp.dI_dE[n]/Ndynamic, ensemble_temp.t);
				
			//costruisco un segnale di cui vado a fare la FFT
			q_p_f[n*Ndynamic + j] = pow( ensemble_temp.q[n]*ensemble_temp.p[n]*(ensemble_temp.dI_dE[n]) , 2.0); // ;
		}
		//ensemble_temp.t += dtdynamic;
	}
	
	double avg_diffusion_temp = 0.0;
	/*
	for(int n=0; n<ensemble_temp.Nparticles; n++)
	{
		//FFT per la traiettoria di ogni particella
		//realft( &q_p_f[n*Ndynamic], Ndynamic );
		for(int j=0; j<Ndynamic; j++)
			avg_diffusion_temp += q_p_f[n*Ndynamic + j];	
	}
	*/
		
	for(int n=0; n<ensemble_temp.Nparticles; n++)
	{
		//FFT per la traiettoria di ogni particella
		realft( &q_p_f[n*Ndynamic], Ndynamic );
		
		avg_diffusion_temp += q_p_f[n*Ndynamic];	
	}
	
	return avg_diffusion_temp /= Ndynamic*ensemble_temp.Nparticles;
	
	//realft( &q_p_f[0], Ndynamic );
	//return q_p_f[0]/Ndynamic;
}

int main(int argc, char *argv[])
{

	int Nensemble = 1000000;
	int nsteps = 10;
	int Ndynamic = 512;
	
	double w = 1.0;
	double k = 1.0;
	double dt = 0.01;
	//double dtdynamic = 0.001;
	double epsilon = 0.2;
	
	double E = 0.1;
	double Emax = w*w*w*w*w*w/(6.0*k*k);
	
	std::cout << "Emax = " << Emax << std::endl;
	std::cout << "Nensemble = " << Nensemble << std::endl;
	std::cout << std::endl;

	double coeff_drift, coeff_diffusion, coeff_diffusion_theoretical;
	std::vector<double> Ncoeff_diffusion, Ntime;
	std::vector<double> Ncoeff_drift;
	
	
	
	Ensemble ensemble(Nensemble, w, k, epsilon, E);
	
	Ensemble ensemble_0 = ensemble;
	
	std::cout << "t iniziale = " << ensemble.t << "\t" << "E avg iniziale = " << ensemble.avg_energy() << std::endl
			  << "azione avg iniziale = " << ensemble.avg_action() << "\t" << "dI_dE[0] iniziale = " << ensemble.dI_dE[0] << std::endl;
	
	double avg_action_0 = ensemble.avg_action();
	
	
	Ensemble ensemble_temp = ensemble;
	double *q_p_f;
	q_p_f = new double[Ndynamic*ensemble_temp.Nparticles];
	
	
	for(int i=0; i<nsteps; i++)
	{
		//Ensemble ensemble_diff = ensemble;
		ensemble.advance(dt);
		//std::cout << "avg action = " << ensemble.avg_action() << std::endl;		
		Ntime.push_back(ensemble.t);
		Ncoeff_drift.push_back(ensemble.avg_action() - avg_action_0);
		Ncoeff_diffusion.push_back( ensemble.avg_action_for_diffusion(ensemble_0) ); 
	}
	
	coeff_drift = slope_linear_regression(Ntime, Ncoeff_drift);//(ensemble.avg_action() - avg_action_0)/dt;
	coeff_diffusion = slope_linear_regression(Ntime, Ncoeff_diffusion);
	coeff_diffusion_theoretical = epsilon*epsilon*theoretical_diffusion(ensemble_temp, Ndynamic, q_p_f);
		
	std::cout << "t finale = " << ensemble.t << "\t" << "E avg finale = " << ensemble.avg_energy() << std::endl
			  << "azione avg finale = " << ensemble.avg_action() << "\t" << "dI_dE[0] finale = " << ensemble.dI_dE[0]<< std::endl;
	
	
	
	std::cout << std::endl << std::endl;
	std::cout << "coeff_drift = " << coeff_drift << std::endl;
	std::cout << "coeff_diffusion = " << coeff_diffusion << std::endl;
	std::cout << "coeff_diffusion_theoretical = " << coeff_diffusion_theoretical << std::endl;
	
	//for(int i=0; i<Ncoeff_diffusion.size(); i++) 		 std::cout << Ncoeff_diffusion[i] 		  << std::endl;
	std::cout << std::endl;
	//for(int i=0; i<Ncoeff_diffusion_dynamic.size(); i++) std::cout << Ncoeff_diffusion_dynamic[i] << std::endl;
	
	int err;
	
	/*
	std::ofstream output;	
	output.open("diffusion_time.txt");
	
	for(int i=0; i<Ncoeff_diffusion.size(); i++)
	{
		output << Ntime[i] << "\t" << Ncoeff_diffusion[i] << std::endl; 	
	}
	
	output.close();
	
	err = system("gnuplot plot_diffusion_time.plt");
	err = system("gnuplot plot_diffusion_time_loglog.plt");
	
	
	std::ofstream spaziofasi;	
	spaziofasi.open("henon_spaziofasi_iniziale.txt");
	
	for(int j=0; j<Ndynamic; j++)
	{		
		for(int n=0; n<ensemble_0.Nparticles; n++)
		{
			ensemble_0.symplectic_advance( ensemble_0.q[n], ensemble_0.p[n], 2.0*PI*ensemble_0.dI_dE[n]/Ndynamic, ensemble_0.t);
			spaziofasi << ensemble_0.q[n] << "\t" << ensemble_0.p[n] << "\t";
		}
		spaziofasi << std::endl;
	}

	spaziofasi.close();
	
	err = system("gnuplot plot_spaziofasi.plt");
	*/
	std::cout << std::endl;
	std::cout << "Fatto!" << std::endl;
	
	return 0;
	
}