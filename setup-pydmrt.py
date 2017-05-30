#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy as np

ext_modules = [ Extension('pydmrt',
                extra_compile_args=['-std=c++11','-fopenmp',"-Wno-sign-compare","-Wno-unused-variable","-Wno-unused-function","-Wno-unused-parameter","-I", "./","-I","../dmrt/dmrt/"],
                extra_link_args=['-lgomp'],
                sources = ['./dmrt/pydmrt.cpp', '../dmrt/dmrt/dmrtreader.cpp', '../dmrt/dmrt/dmrtalg2.cpp',
                           '../dmrt/dmrt/dmrtmain.cpp']) ]

setup(
        name = 'pydmrt',
        version = '0.1',
        include_dirs = [np.get_include()], #Add Include path of numpy
        ext_modules = ext_modules
      )
