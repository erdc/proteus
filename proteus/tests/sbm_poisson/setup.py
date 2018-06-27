import sys
import setuptools
from distutils.core import setup, Extension

import numpy
from Cython.Distutils import build_ext

## \file setup.py setup.py
#  \brief The python script for building proteus
#
#  Set the DISTUTILS_DEBUG environment variable to print detailed information while setup.py is running.
#

from proteus import config
from proteus.config import *

###to turn on debugging in c++
##\todo Finishing cleaning up setup.py/setup.cfg, config.py...
from distutils import sysconfig
cv = sysconfig.get_config_vars()
cv["OPT"] = cv["OPT"].replace("-DNDEBUG","-DDEBUG")
cv["OPT"] = cv["OPT"].replace("-O3","-g")
cv["CFLAGS"] = cv["CFLAGS"].replace("-DNDEBUG","-DDEBUG")
cv["CFLAGS"] = cv["CFLAGS"].replace("-O3","-g")
cv["CFLAGS"]+=" -march=native"


PROTEUS_PETSC_EXTRA_LINK_ARGS = getattr(config, 'PROTEUS_PETSC_EXTRA_LINK_ARGS', [])
PROTEUS_PETSC_EXTRA_COMPILE_ARGS = getattr(config, 'PROTEUS_PETSC_EXTRA_COMPILE_ARGS', [])

proteus_install_path = os.path.join(sysconfig.get_python_lib(), 'proteus')

# handle non-system installations
for arg in sys.argv:
    if arg.startswith('--root'):
        proteus_install_path = proteus_install_path.partition(sys.prefix + '/')[-1]
        break
    if arg.startswith('--prefix'):
        proteus_install_path = proteus_install_path.partition(sys.prefix + '/')[-1]
        break

setup(name='ChMBDModel',
      version='0.0.1',
      description='Python tools for multiphysics modeling',
      author='Chris Kees, Matthew Farthing, et al.',
      author_email='chris.kees@us.army.mil',
      url='http://proteus.usace.army.mil',
      cmdclass = {'build_ext':build_ext},
      ext_modules=[
#                     Extension("Chrono",['Chrono.pyx'],
#                              depends=['Chrono.h'],
#                              language='c++',
#                              include_dirs=['/home/jovyan/proteus/linux2/include/chrono/collision/convexdecomposition/HACD',
#                                            '/home/jovyan/proteus/linux2/include/chrono/collision/bullet',
#                                            numpy.get_include(),'proteus',config.PROTEUS_INCLUDE_DIR],
#                              library_dirs=[config.PROTEUS_LIB_DIR,config.PROTEUS_LIB_DIR[:-3]+'lib64'],
#                              libraries=['ChronoEngine',
#                                         'stdc++','m'],
#                              extra_compile_args=["-std=c++11"]),
                    Extension("cPoisson_M1",["cPoisson_M1.pyx"],
                              depends=["Poisson_M1.h"] + ["ModelFactory.h","CompKernel.h"],
                              language="c++",
                              extra_compile_args=PROTEUS_OPT,
                              include_dirs=[numpy.get_include(),
                                            '/Users/yy/p/proteus/proteus'],
                              library_dirs=[config.PROTEUS_LIB_DIR,config.PROTEUS_LIB_DIR[:-3]+'lib64'],
                                ),
                ]
      )
