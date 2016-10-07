#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <stdlib.h>

#include "cranknicolson.h"

#define PI 3.14159265359
#define TWOPI	(2.0*PI)

//gaussian distribution
double gauss(double x, double mean, double devstd)
{
	return exp( -0.5*pow( (x - mean)/devstd , 2.0))/(devstd*sqrt(TWOPI));
}

int main(void)
{
	//store values loaded from file
	std::vector<double> numbers;

	std::ifstream in("diffusion_action.txt"); //Create an input file stream

	double number;  //Variable to hold each number as it is read
	
    while (in >> number) { //Read number using the extraction (>>) operator
		numbers.push_back(number);//Add the number to the end of the array
	}

	in.close(); //Close the file stream

	
	int N = numbers.size()/2 - 2; //without boundaries that are used to evaluate D in the edges
			
	//memory allocation
	double *x,*ro,*D,*mu;
	x  = new double[N+2]; 
	ro = new double[N];
	D  = new double[N+2];
	mu = new double[N+2];
	
	//initialize x, ro, D, mu 
	for(int i=0; i<N+2; i++)
	{
		x[i] = numbers[2*i];
		D[i] = numbers[2*i+1];
		//mu[i] = -0.5*x[i]; //not used
	}
	
	//ro initialized as gaussian distribution
	for(int i=1; i<N; i++)
	{
		ro[i] = gauss(x[i], 0.040, 0.005);
	} 
	//boundary conditions
	ro[0] = 0.0;
	ro[N-1] = 0.0;

	int Npassi = 20000;
	double dt = 0.0001;
	double h = (x[N+1]-x[0])/(double)(N+1);
	std::cout << "h = " << h << std::endl;

	std::cout << std::endl;
	std::cout << "Maximum time step = " << h*h/D[N+1] << std::endl;
		
	
	//for(int i=0; i<N; i++) std::cout << D[i] << "\t" << x[i] << "\t" << ro[i] << std::endl;
	
	//compute numerically the integral of the initial distribution  
	double sum = 0.0;
	for(int i=0; i<N; i++) sum += h*ro[i];
	std::cout << "sum iniziale = " << sum << std::endl;
	
	
	//apply crank nicolson and find the final distribution resulted from the diffusion process
	double *res = cranknicolson2(D, ro, N, Npassi, h, dt);
	//double *res = cranknicolson_fokkerplanck(mu, D, ro, N, Npassi, h, dt);
	
	//for(int i=0; i<N; i++) std::cout << x[i] << "\t" << ro[i] << "\t" << res[i] << std::endl;
	
	//compute numerically the integral of the final distribution  
	sum = 0.0;
	for(int i=0; i<N; i++) sum += h*res[i];
	std::cout << "sum finale = " << sum << std::endl;
	
	
	//save to file and plot
	int err;
	
	std::ofstream output;	
	output.open("henon_diffusion_cranknicolson.txt");
	
	for(int i=0; i<N; i++)
	{
		output << x[i] << "\t" << ro[i] << "\t" << res[i] << std::endl; 	
	}
	
	output.close();
	
	err = system("gnuplot plot_henon_diffusion_cranknicolson.plt");
	
	return 0;	
}