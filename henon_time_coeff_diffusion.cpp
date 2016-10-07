#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm> // std::max
#include <math.h>	// sqrt
#include <complex> //std::sqrt

#include "ensemble_henon_smart.h"

#define PI 3.14159265359

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

double theoretical_diffusion(Ensemble ensemble_temp, int Ndynamic)
{	
	double *q_p_f;
	q_p_f = new double[Ndynamic*ensemble_temp.Nparticles];
	
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
	//valore medio senza fft
	for(int n=0; n<ensemble_temp.Nparticles; n++)
		for(int j=0; j<Ndynamic; j++)
			avg_diffusion_temp += q_p_f[n*Ndynamic + j];	
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

	int Nensemble = 10000;
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
	
	//save copy of the initial state of the ensemble
	Ensemble ensemble_0 = ensemble;
	//output average parameters of the ensemble
	std::cout << "t iniziale = " << ensemble.t << "\t" << "E avg iniziale = " << ensemble.avg_energy() << std::endl
			  << "azione avg iniziale = " << ensemble.avg_action() << "\t" << "dI_dE avg iniziale = " << ensemble.avg_dI_dE() << std::endl;
	
	//initial average  action of the ensemble
	double avg_action_0 = ensemble.avg_action();
	//temporary copy of the initial ensemble for the theoretical diffusion
	Ensemble ensemble_temp = ensemble;

	for(int i=0; i<nsteps; i++)
	{
		ensemble.advance(dt);
		
		Ntime.push_back(ensemble.t);
		Ncoeff_drift.push_back(ensemble.avg_action() - avg_action_0);
		Ncoeff_diffusion.push_back( ensemble.avg_action_for_diffusion(ensemble_0) ); 
	}
	
	coeff_drift 	= slope_linear_regression(Ntime, Ncoeff_drift);
	coeff_diffusion = slope_linear_regression(Ntime, Ncoeff_diffusion);
	coeff_diffusion_theoretical = epsilon*epsilon*theoretical_diffusion(ensemble_temp, Ndynamic);
		
	std::cout << "t finale = " << ensemble.t << "\t" << "E avg finale = " << ensemble.avg_energy() << std::endl
			  << "azione avg finale = " << ensemble.avg_action() << "\t" << "dI_dE avg iniziale = " << ensemble.avg_dI_dE() << std::endl;
	
	
	
	std::cout << std::endl << std::endl;
	std::cout << "coeff_drift = " << coeff_drift << std::endl;
	std::cout << "coeff_diffusion = " << coeff_diffusion << std::endl;
	std::cout << "coeff_diffusion_theoretical = " << coeff_diffusion_theoretical << std::endl;
	
	//for(int i=0; i<Ncoeff_diffusion.size(); i++) 		 std::cout << Ncoeff_diffusion[i] 		  << std::endl;
	std::cout << std::endl;
	//for(int i=0; i<Ncoeff_diffusion_dynamic.size(); i++) std::cout << Ncoeff_diffusion_dynamic[i] << std::endl;
	
	int err;
	
	std::ofstream output;	
	output.open("diffusion_time.txt");
	
	for(int i=0; i<Ncoeff_diffusion.size(); i++)
	{
		output << Ntime[i] << "\t" << Ncoeff_diffusion[i] << std::endl; 	
	}
	
	output.close();
	
	err = system("gnuplot plot_diffusion_time.plt");
	//err = system("gnuplot plot_diffusion_time_loglog.plt");
	
	std::cout << std::endl;
	std::cout << "Fatto!" << std::endl;
	
	return 0;
	
}