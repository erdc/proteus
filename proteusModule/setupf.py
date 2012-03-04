from distutils.core import setup, Extension
import numpy

from distutils import sysconfig
from numpy.distutils.core import setup, Extension
## \file setup.py setup.py
#  \brief The python script for building proteus
#
#  Set the DISTUTILS_DEBUG environment variable to print detailed information while setup.py is running.
#
try:
    from config import *
except:
    raise RuntimeError("Missing or broken config.py file. See proteusConfig for examples")

cv = sysconfig.get_config_vars()
cv["OPT"] = cv["OPT"].replace("-DNDEBUG","-DDEBUG")
cv["OPT"] = cv["OPT"].replace("-O3","")
cv["CFLAGS"] = cv["CFLAGS"].replace("-DNDEBUG","-DDEBUG")
cv["CFLAGS"] = cv["CFLAGS"].replace("-O3","")

setup(name='ftracking',
      version='0.0.1',
      description='Fortran based Python tools for particle tracking',
      author='Pearce Cheng and Matthew Farthing',
      author_email='matthew.w.farthing@usace.army.mil',
      package_dir={'proteus':'src'},
      ext_package='proteus',
      ext_modules=[Extension('ftracking',
                             ['src/ftracking.f'],
                             include_dirs=[numpy.get_include(),'include'],
                             libraries=['m'],
                             extra_link_args=['-g']+PROTEUS_EXTRA_FC_LINK_ARGS,
                             extra_compile_args=['-g']+PROTEUS_EXTRA_FC_COMPILE_ARGS)],
      requires=['numpy']
      )
