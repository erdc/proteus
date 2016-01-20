import numpy
from distutils.core import setup
from distutils.extension import Extension

USE_CYTHON= True

ext = '.pyx' if USE_CYTHON else '.c'

extensions_pyx = [Extension("cev_utils",['cev_utils'+ext],
                        include_dirs=[numpy.get_include()])]

ext_extensions = []
if USE_CYTHON:
    from Cython.Build import cythonize
    ext_extensions += cythonize(extensions_pyx) 
setup(name='ev_utils',
      description='Utilities for Entropy Viscosity approximations',
      ext_modules=ext_extensions,
      requires=['numpy','scipy']
      )

