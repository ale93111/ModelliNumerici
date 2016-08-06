# -*- coding: utf-8 -*-
"""
Created on Sat Aug  6 18:21:37 2016

@author: alessandro
"""

import ensemble 
import numpy as np
import matplotlib.pyplot as plt
import time
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
ef = ensemble.PyEnsemble(Nensemble,w,k,epsilon,mean,devstd,todo)

#%%
t = time.time()
for i in range(nsteps):
    ef.pyadvance(dt)
    
print(time.time() - t)
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
plt.figure(figsize=(12,12))
plt.hist(ei.q, bins=30, color='blue', normed=1, alpha=0.5)
plt.hist(ef.q, bins=30, color='red', normed=1, alpha=0.5)
plt.xlabel('q')
plt.ylabel('# Occurrences')
#%%
plt.figure(figsize=(12,12))
plt.hist(ei.p, bins=30, color='blue', normed=1, alpha=0.5)
plt.hist(ef.p, bins=30, color='red', normed=1, alpha=0.5)
plt.xlabel('p')
plt.ylabel('# Occurrences')
#%%
plt.figure(figsize=(12,12))
plt.hist(ei.E, bins=30, color='blue', normed=1, alpha=0.5)
plt.hist(ef.E, bins=30, color='red', normed=1, alpha=0.5)
plt.xlabel('Energy')
plt.ylabel('# Occurrences')
#%%
plt.figure(figsize=(12,12))
plt.hist(ei.action, bins=30, color='blue', normed=1, alpha=0.5)
plt.hist(ef.action, bins=30, color='red', normed=1, alpha=0.5)
plt.xlabel('Action')
plt.ylabel('# Occurrences')
#%%
plt.figure(figsize=(12,12))
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
        
#%%
temparray = np.array(arraylist)
array = np.reshape(temparray, (3, int(np.size(temparray)/3)))

for i in range(int(np.size(array)/3)):
    array[0][i] = arraylist[i][0] #x
    array[1][i] = arraylist[i][1] #y initial
    array[2][i] = arraylist[i][2] #y final
#%%
plt.figure(figsize=(12,12))
plt.plot(array[0][:], array[1], 'b--',linewidth=4.0, label='Initial distribution - Cranck-nicolson')
plt.plot(array[0][:], array[2], 'r--',linewidth=4.0, label='Final distribution - Cranck-nicolson')
ni, binsi, patchesi = plt.hist(ei.action, bins=30, color='blue', normed=1, alpha=0.5, label='Initial distribution - MonteCarlo')
nf, binsf, patchesf = plt.hist(ef.action, bins=30, color='red', normed=1, alpha=0.5, label='Final distribution - MonteCarlo')

plt.xlabel('Action')
plt.ylabel('Normalized occurrences')
plt.legend()

