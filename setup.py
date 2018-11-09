#!/usr/bin/env python

from __future__ import print_function

import sys, os, platform

from setuptools import setup
from setuptools import Extension
import distutils.sysconfig

from Cython.Build import cythonize

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

      packages = ['AAKwrapper'],
      package_dir = {'AAKwrapper': 'AAKwrapper'},

      ext_modules = cythonize(Extension('AAKwrapper.AAKwrapper',['AAKwrapper/AAKwrapper.pyx'],
                                        language = "c++",
                                        include_dirs = ['./include', numpy.get_include()],
                                        libraries = ['gslcblas', 'gsl', 'KS', 'IEKG', 'LB', 'NR', 'RRGW', 'GKG', 'Circ'],
                                        library_dirs = ['/usr/local/lib', './lib'],
                                        extra_compile_args = ["-Wno-unused-function"],
                                        extra_link_args = linkArgs))
                              )

