# ModelliNumerici

All scripts need to be compiled with the "-std=c++11" flag, also the optional flag "-fopenmp" to use OpenMP

The diffusion coefficients values are computed using "henon_find_coeff_diffusion.cpp". 
These values are used by "henon_diffusion_cranknicolson.cpp".
Finally "henon_diffusion_cranknicolson.cpp" and "henon_noise_diffusion.cpp" results can be compared.