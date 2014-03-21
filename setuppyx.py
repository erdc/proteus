from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy,os
try:
    from config import *
except:
    raise RuntimeError("Missing or broken config.py file. See proteusConfig for examples")

setup(
    cmdclass = {'build_ext':build_ext},
    ext_package='proteus',
    package_dir={'proteus':'src'},
    ext_modules = [Extension("waveFunctions",['src/waveFunctions.pyx','src/transportCoefficients.c'],
                            include_dirs=[numpy.get_include(),'include']),
                   Extension("subsurfaceTransportFunctions",['src/subsurfaceTransportFunctions.pyx'],
                            include_dirs=[numpy.get_include(),'include'])]
    )
