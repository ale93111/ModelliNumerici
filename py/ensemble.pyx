# -*- coding: utf-8 -*-
"""
Created on Sat Aug  6 23:31:53 2016

@author: alessandro
"""
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref

cdef extern from "/home/alessandro/ModelliNumerici/ensemble_henon_smart.h":
    
    void find_roots(double & root1, double & root2, double & root3, const double a, const double b, const double c)
    
    cdef cppclass Ensemble:
        
        vector[double] p, q, E, action, dI_dE
        vector[int] deletedlist
        double w, k, epsilon, t
        int Nparticles
	
 
        double avg_energy()
        double avg_action()
        double avg_dI_dE()
        double energy_max()
        void check_energy(int index)
        double energy(const double & q, const double & p)
        double pex(const double & E, const double & x)
        double tex(const double & E, const double & x) 
        double action_elliptic(double root1, double root2, double root3)
        double dI_dE_elliptic(double root1, double root2, double root3)
        void symplectic_advance(double & q, double & p, const double dt, const double t)
        void advance( const double dt )
        void newtonfunc(double x, double *f, double *df, double I0)
        double rtsafe(double x1, double x2, double xacc, double offset)
        void allocator( int N )
        
        Ensemble() except +
        Ensemble(int N, double wi, double ki, double epsiloni, double Ei) except +
        Ensemble(int N, double wi, double ki, double epsiloni, double actioni, char Vi) except +
        Ensemble(int N, double wi, double ki, double epsiloni, double mean, double devstd, char Vi) except +
        Ensemble& operator=(Ensemble&& other) except +
        Ensemble& operator=(Ensemble& other) except + 
        

    
cdef class PyEnsemble:
    cdef Ensemble c_ensemble      # hold a C++ instance which we're wrapping
    
    def __cinit__(self, int N=0, double w=0, double k=0, double epsilon=0, double EImean=0, double devstd=0, char Vi='0'):
        if N == 0:
            self.c_ensemble = Ensemble()
        if devstd == 0 and Vi == '0':
            self.c_ensemble = Ensemble(N,w,k,epsilon,EImean)
        if devstd == 0 and Vi != '0':
            self.c_ensemble = Ensemble(N,w,k,epsilon,EImean,Vi)
        if devstd != 0 and Vi != '0':
            self.c_ensemble = Ensemble(N,w,k,epsilon,EImean,devstd,Vi)
        
    def pyavg_energy(self):
        return self.c_ensemble.avg_energy()
    def pyavg_action(self):
        return self.c_ensemble.avg_action()
    def pyavg_dI_dE(self):
        return self.c_ensemble.avg_dI_dE()
    def pyenergy_max(self):
        return self.c_ensemble.energy_max()
        
    def pyadvance(self, dt):
        self.c_ensemble.advance(dt)
    

    @property
    def Nparticles(self):
        return self.c_ensemble.Nparticles    
    @property
    def p(self):
        return self.c_ensemble.p
    @property
    def q(self):
        return self.c_ensemble.q
    @property
    def E(self):
        return self.c_ensemble.E
    @property
    def action(self):
        return self.c_ensemble.action
    @property
    def dI_dE(self):
        return self.c_ensemble.dI_dE
    @property
    def deletedlist(self):
        return self.c_ensemble.deletedlist
    @property
    def t(self):
        return self.c_ensemble.t
