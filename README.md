# ModelliNumerici

# Introduction
This work was done to study the diffusion in the dynamics of particle beams in accelerator
physics. The phenomenon can be described by the Hénon map whose linear frequency is
stochastically perturbed. Finally, a theoretical model for the action diffusion based on the
Fokker-Planck equation is introduced and its solution agrees with the simulation process.


Read report.pdf for a more detailed explanation.

All scripts need to be compiled with the "-std=c++11" flag, also the optional flag "-fopenmp" to use OpenMP and/or "-O3" for optimization.
The main part of the code was also integrated in Python 3+ using Cython. A python script "ensemble_setup.py" is available to install the package.

The diffusion coefficients values are computed using "henon_find_coeff_diffusion.cpp". 
These values are used by "henon_diffusion_cranknicolson.cpp".
Finally "henon_diffusion_cranknicolson.cpp" and "henon_noise_diffusion.cpp" results can be compared.



