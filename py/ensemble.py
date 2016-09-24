# -*- coding: utf-8 -*-
"""
Created on Sat Aug  6 18:21:37 2016

@author: alessandro
"""

import ensemble 
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy as sp

sizefig = 7
#%%
Nensemble = 10000
w = 1.0
k = 1.0
epsilon = 0.1;
mean = 0.040
devstd = 0.005

nsteps = 2000
dt = 0.01

todo = 1
#%%
ei = ensemble.PyEnsemble(Nensemble,w,k,epsilon,mean,devstd,todo)
#temporary workaround to a python copy constructor
ef = ensemble.PyEnsemble(0,0,0,0,0,0,0,ei)


#%%
#t = time.time()
%%time
for i in range(nsteps):
    ef.pyadvance(dt)
    
#print("time=",time.time() - t)
#%%
print("ENSEMBLE INIZIALE")
print(" N particles = ", ei.Nparticles    ,'\n',
      "time = "        , ei.t             ,'\n',    
      "max energy = "  , ei.pyenergy_max(),'\n','\n', 
      "avg energy = "  , ei.pyavg_energy(),'\n',
      "avg action = "  , ei.pyavg_action(),'\n',
      "avg dI_dE = "   , ei.pyavg_dI_dE() ,'\n')
print("ENSEMBLE FINALE")
print(" N particles = ", ef.Nparticles    ,'\n',
      "time = "        , ef.t             ,'\n',    
      "max energy = "  , ef.pyenergy_max(),'\n','\n', 
      "avg energy = "  , ef.pyavg_energy(),'\n',
      "avg action = "  , ef.pyavg_action(),'\n',
      "avg dI_dE = "   , ef.pyavg_dI_dE() ,'\n')
      
#%%
plt.figure(figsize=(sizefig,sizefig))
plt.hist(ei.q, bins=30, color='blue', normed=1, alpha=0.5)
plt.hist(ef.q, bins=30, color='red', normed=1, alpha=0.5)
plt.xlabel('q')
plt.ylabel('# Occurrences')
#%%
plt.figure(figsize=(sizefig,sizefig))
plt.hist(ei.p, bins=30, color='blue', normed=1, alpha=0.5)
plt.hist(ef.p, bins=30, color='red', normed=1, alpha=0.5)
plt.xlabel('p')
plt.ylabel('# Occurrences')
#%%
plt.figure(figsize=(sizefig,sizefig))
plt.hist(ei.E, bins=30, color='blue', normed=1, alpha=0.5)
plt.hist(ef.E, bins=30, color='red', normed=1, alpha=0.5)
plt.xlabel('Energy')
plt.ylabel('# Occurrences')
#%%
plt.figure(figsize=(sizefig,sizefig))
plt.hist(ei.action, bins=30, color='blue', normed=1, alpha=0.5)
plt.hist(ef.action, bins=30, color='red', normed=1, alpha=0.5)
plt.xlabel('Action')
plt.ylabel('# Occurrences')
#%%
plt.figure(figsize=(sizefig,sizefig))
plt.hist(ei.dI_dE, bins=30, color='blue', normed=1, alpha=0.5)
plt.hist(ef.dI_dE, bins=30, color='red', normed=1, alpha=0.5)
plt.xlabel('dI_dE')
plt.ylabel('# Occurrences')
#%%
path = '/home/alessandro/ModelliNumerici/'
fname = 'henon_diffusion_cranknicolson.txt'
with open(path+fname) as f:
    #w, h, z = [float(x) for x in next(f).split()] # read first line
    arraylist = []
    for line in f: # read rest of lines
        arraylist.append([float(x) for x in line.split()])
        

temparray = np.array(arraylist)
array = np.reshape(temparray, (3, int(np.size(temparray)/3)))

for i in range(int(np.size(array)/3)):
    array[0][i] = arraylist[i][0] #x
    array[1][i] = arraylist[i][1] #y initial
    array[2][i] = arraylist[i][2] #y final
#%%
plt.figure(figsize=(sizefig,sizefig))
plt.plot(array[0][:], array[1], 'b--',linewidth=4.0, label='Initial distribution - Cranck-nicolson')
plt.plot(array[0][:], array[2], 'r--',linewidth=4.0, label='Final distribution - Cranck-nicolson')
ni, binsi, patchesi = plt.hist(ei.action, bins=30, color='blue', normed=1, alpha=0.5, label='Initial distribution - MonteCarlo')
nf, binsf, patchesf = plt.hist(ef.action, bins=100, color='red', normed=1, alpha=0.5, label='Final distribution - MonteCarlo')

plt.xlabel('Action')
plt.ylabel('Normalized occurrences %')
plt.legend()
#%%

Emax = 0.16666666
Nensemble = 9
energy = sp.linspace(0.1*Emax, 0.9*Emax, Nensemble)

#TEST SPAZIO FASI
test = []
for i in range(Nensemble):
    test.append(ensemble.PyEnsemble(1,w,k,epsilon,energy[i]))

#%%
npoints = 512
dt = []
for i in range(Nensemble):
    dt.append( 2.1*np.pi*test[i].dI_dE[0]/npoints)


p = []
q = []

for i in range(Nensemble):
    
    test_p = []
    test_q = []
    for j in range(npoints):
        test[i].pysymplectic_advance(0, dt[i])
        test_q.append(test[i].q[0])
        test_p.append(test[i].p[0])
    
    p.append(test_p)
    q.append(test_q)

    
#%%
plt.figure(figsize=(sizefig,sizefig))
for i in range(Nensemble):
    plt.plot(q[i], p[i], 'b', linewidth=1.0)
plt.plot(0,0,'.')

plt.axis([-0.6, 1.0, -0.8, 0.8])
plt.xlabel('q')
plt.ylabel('p')
#%%
#phase space distribution
#TODO (low priority) non è corretta, dovrebbe essere distribuita uniforme rispetto l'angolo
Nparticles = 100
test = ensemble.PyEnsemble(Nparticles,w,k,epsilon,0.1)

plt.figure(figsize=(sizefig,sizefig))
for i in range(Nparticles):
    plt.plot(test.q, test.p, '.')

plt.axis([-0.6, 1.0, -0.8, 0.8])
plt.xlabel('q')
plt.ylabel('p')
