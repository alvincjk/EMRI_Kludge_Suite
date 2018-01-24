#!/usr/bin/env python

from __future__ import print_function

import sys, os, platform

from setuptools import setup
from setuptools import Extension
import distutils.sysconfig

import numpy

# maybe need rpath links to shared libraries on Linux
# if platform.system() == 'Linux':
#    linkArgs = ['-Wl,-R{}/lib'.format(...)]
#else:
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

      ext_modules = [Extension('AAKwrapper/_AAKwrapper',['AAKwrapper/AAKwrapper.i'],
                               language = 'c++',
                               swig_opts = ['-c++'],
                               include_dirs = ['./include', numpy.get_include()],
                               libraries = ['gslcblas', 'gsl', 'KS', 'IEKG', 'LB', 'NR', 'RRGW', 'GKG', 'Circ'],
                               library_dirs = ['/usr/local/lib', './lib'],
                               extra_compile_args = ["-Wno-unused-function"],
                               extra_link_args = linkArgs)]
      )
