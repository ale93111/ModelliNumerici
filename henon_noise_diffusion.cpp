#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm> // std::max
#include <math.h>	// sqrt
#include <complex> //std::sqrt

//#include"ensemble_henon.h"
#include"ensemble_henon_smart.h"

#define PI 3.14159265359

int main(int argc, char *argv[])
{

	int Nensemble = 10000;
	int nsteps = 10;
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
	/*
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
	*/
	return 0;
	
}