#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm> // std::max
#include <math.h>	// sqrt
#include <complex> //std::sqrt

#include "ensemble_henon_smart.h"
#include "cranknicolson.h"

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

	for(int n=0; n<ensemble_temp.Nparticles; n++)
	{
		//FFT per la traiettoria di ogni particella
		realft( &q_p_f[n*Ndynamic], Ndynamic );
		
		avg_diffusion_temp += q_p_f[n*Ndynamic];	
	}
	
	return avg_diffusion_temp /= Ndynamic*ensemble_temp.Nparticles;
}

int main(void)
{
	double w = 1.0;
	double k = 1.0;
	int Nensemble = 100;
	int Ndynamic = 512;
	double epsilon = 0.2;
	
	double Emax = w*w*w*w*w*w/(6.0*k*k);
	
	std::cout << "Emax = " << Emax << std::endl;
	std::cout << "Nensemble = " << Nensemble << std::endl;
	std::cout << std::endl;
	
	int n = 201; //add +1 for the last element of the interval
		
	double Ea, Eb;
	
	Ea = 0.1*Emax;
	Eb = 0.9*Emax;
	
	Ensemble *ensemble = new Ensemble[n];
	
	ensemble[0]   = Ensemble(Nensemble, w, k, epsilon, Ea);
	ensemble[n-1] = Ensemble(Nensemble, w, k, epsilon, Eb);
	
	double *I = new double[n];
	double Ia = ensemble[0].avg_action();
	double Ib = ensemble[n-1].avg_action();
	
	double dI = (Ib - Ia)/(double)(n-1);
	std::cout << "dI = " << dI << std::endl;
	I[0] = Ia;
	I[n-1] = Ib;
	for(int i=1; i<n; i++) I[i] = I[i-1] + dI;
	
	//for(int i=0; i<n; i++) std::cout << "I [" << i << "] = " << I[i] << std::endl;
	
	for(int i=1; i<n-1; i++) ensemble[i] = Ensemble(Nensemble, w, k, epsilon, I[i], 'i');

	double *Ncoeff_diffusion = new double[n];
	
	for(int i=0; i<n; i++) Ncoeff_diffusion[i] = epsilon*epsilon*theoretical_diffusion(ensemble[i], Ndynamic);
	
	//for(int i=0; i<n; i++) std::cout << "Coeff_diffusion [" << i << "] = " << Ncoeff_diffusion[i] << std::endl;
	
	
	//save to file and plot
	int err;
	std::ofstream output;	
	output.open("diffusion_action.txt");
	
	for(int i=0; i<n; i++)
	{
		output << I[i] << "\t" << Ncoeff_diffusion[i] << std::endl; 	
	}
	
	output.close();
	
	err = system("gnuplot plot_diffusion_action.plt");
	
	
	
	return 0;
}