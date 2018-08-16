#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy as np

from distutils.command.build_ext import build_ext
from distutils.sysconfig import customize_compiler


class my_build_ext(build_ext):
    def build_extensions(self):
        customize_compiler(self.compiler)
        try:
            self.compiler.compiler_so.remove("-Wstrict-prototypes")
        except (AttributeError, ValueError):
            pass
        build_ext.build_extensions(self)

ext_modules = [ Extension('pydmrt',
                extra_compile_args=['-std=c++11',"-Wno-sign-compare","-Wno-unused-function","-Wno-unused-parameter","-I", "./","-I","../dmrt/dmrt/"],
                sources = ['./dmrt/pydmrt.cpp', '../dmrt/dmrt/dmrtreader.cpp', '../dmrt/dmrt/dmrtalg2.cpp',
                           '../dmrt/dmrt/dmrtmain.cpp']) ]

setup(
        cmdclass = {'build_ext': my_build_ext},
        name = 'pydmrt',
        version = '0.1',
        include_dirs = [np.get_include()], #Add Include path of numpy
        ext_modules = ext_modules
      )
