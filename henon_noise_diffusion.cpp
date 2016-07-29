#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm> // std::max
#include <math.h>	// sqrt
#include <complex> //std::sqrt

//#include"fourier.h"
//#include"elliptic_int.h"
//#include"ensemble_henon.h"
#include"ensemble_henon_smart.h"

#define npoints 1000
#define niterate 10
#define PI 3.14159265359


double theoretical_diffusion(Ensemble ensemble_temp, int Ndynamic/*, double *q_p_f*/)
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
}

int main(int argc, char *argv[])
{

	int Nensemble = 10000;
	int nsteps = 2000;
	int Ndynamic = 512;
	
	double w = 1.0;
	double k = 1.0;
	//double dtdynamic = 0.001;
	double epsilon = 0.1;
	double dt = 0.01;
	
	double mean = 0.040;
	double devstd = 0.005;
	double E = 0.1;
	double Emax = w*w*w*w*w*w/(6.0*k*k);
	
	std::cout << "Emax = " << Emax << std::endl;
	std::cout << "Nensemble = " << Nensemble << std::endl;
	std::cout << std::endl;

	//Ensemble ensemble(Nensemble, w, k, epsilon, E);
	Ensemble ensemble(Nensemble, w, k, epsilon, mean, devstd, 'i');
	
	Ensemble ensemble_0 = ensemble; //save a copy of the initial ensemble
	
	std::cout << "t iniziale = " << ensemble.t << "\t" << "E avg iniziale = " << ensemble.avg_energy() << std::endl
			  << "azione avg iniziale = " << ensemble.avg_action() << "\t" << "dI_dE avg iniziale = " << ensemble.avg_dI_dE() << std::endl;
	
	
	for(int i=0; i<nsteps; i++)
	{
		ensemble.advance(dt);
	}
	
	std::cout << "t finale = " << ensemble.t << "\t" << "E avg finale = " << ensemble.avg_energy() << std::endl
			  << "azione avg finale = " << ensemble.avg_action() << "\t" << "dI_dE avg finale = " << ensemble.avg_dI_dE()<< std::endl;
	
	
	std::cout << std::endl;
	
	int err;
	
	//salvare su file
	std::ofstream output;	
	output.open("henon_noise_diffusion.txt");
	
	for(int i=0; i<ensemble.Nparticles; i++)
	{
		output << ensemble_0.action[i] << "\t" << ensemble.action[i]  << std::endl; 	
	}
	
	output.close();
	
	err = system("gnuplot plot_henon_noise_diffusion.plt");

	std::cout << std::endl;
	std::cout << "Fatto!" << std::endl;
	
	return 0;
	
}