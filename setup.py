#!/usr/bin/env python

from __future__ import print_function

import sys, os, platform

from setuptools import setup
from setuptools import Extension
import distutils.sysconfig

from Cython.Build import cythonize

import numpy

# maybe need rpath links to shared libraries on Linux
if platform.system() == 'Linux':
    linkArgs = ['-Wl,-R{}/lib'.format(tempo2)]
else:
    linkArgs = []

setup(name = 'AAKwrapper',
      version = '0.1',
      description = 'A Python wrapper for AAK',

      # author = 'Michele Vallisneri',
      # author_email = 'vallis@vallis.org',
      # url = 'https://github.com/ajchua/libstempo',

      packages = ['AAKwrapper'],
      package_dir = {'AAKwrapper': 'AAKwrapper'},
      # package_data = {'libstempo': ['data/*', 'ecc_vs_nharm.txt']},

      # py_modules = ['libstempo.like','libstempo.multinest','libstempo.emcee'],

      ext_modules = cythonize(Extension('AAKwrapper.AAKwrapper',['AAKwrapper/AAKwrapper.pyx'],
                                        language = "c++",
                                        include_dirs = ['./include', numpy.get_include()],
                                        libraries = ['gslcblas', 'gsl', 'Circ', 'GKG', 'IEKG', 'KS', 'LB', 'NR', 'RRGW'],
                                        library_dirs = ['/usr/local/lib', './lib'],
                                        extra_compile_args = ["-Wno-unused-function"],
                                        extra_link_args = linkArgs))
                              )
