make:
	g++ henon_time_coeff_diffusion.cpp    -std=c++11 -fopenmp -O3 -o time_coeff_diffusion.out
	g++ henon_noise_diffusion.cpp         -std=c++11 -fopenmp -O3 -o noise_diffusion.out
	g++ henon_diffusion_cranknicolson.cpp -std=c++11 -o diffusion_cranknicolson.out
	g++ henon_find_coeff_diffusion.cpp    -std=c++11 -o find_coeff_diffusion.out


