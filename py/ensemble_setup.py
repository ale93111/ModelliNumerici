# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 01:08:45 2016

@author: alessandro
"""
#from distutils.core import setup
#from Cython.Build import cythonize
#
#setup(ext_modules = cythonize(
#           "/home/alessandro/ModelliNumerici/py/ensemble.pyx",                 # our Cython source
#           #sources=["/home/alessandro/ModelliNumerici/fourier.h,/home/alessandro/ModelliNumerici/elliptic_int.h"],  # additional source file(s)
#           language="c++",             # generate C++ code
#      ))

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    name = 'Ensemble',
    ext_modules=[
    Extension('ensemble',
              sources=['ensemble.pyx'],
              extra_compile_args=['-std=c++11','-O3','-fopenmp'],
              extra_link_args=['-lgomp'],
              language='c++')
    ],
    cmdclass = {'build_ext': build_ext}
)

#%%
