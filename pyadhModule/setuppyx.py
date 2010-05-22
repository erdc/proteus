from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy,os
try:
    from config import *
except:
    print "Using defaultConfig.py"
    from defaultConfig import *
setup(
    cmdclass = {'build_ext':build_ext},
    ext_package='pyadh',
    ext_modules = [Extension("waveFunctions",['pyadh/waveFunctions.pyx','pyadh/transportCoefficients.c'],
                            include_dirs=[numpy.get_include(),'include']),
                   Extension("subsurfaceTransportFunctions",['pyadh/subsurfaceTransportFunctions.pyx'],
                            include_dirs=[numpy.get_include(),'include']),
                   Extension("pskRelations",['pyadh/pskRelations.pyx'],
                             include_dirs=[numpy.get_include(),'include']),
                   Extension('cadh',
                             ['pyadh/cadh.pyx','pyadh/cadhimpl.c'],
                             define_macros=[('PYADH_SUPERLU_H',PYADH_SUPERLU_H),
                                            ('PYADH_LAPACK_H',PYADH_LAPACK_H),
                                            ('PYADH_LAPACK_INTEGER',PYADH_LAPACK_INTEGER),
                                            ('PYADH_BLAS_H',PYADH_BLAS_H),
                                            ('_DEBUG',1),
                                            ('_MESSG',1),
                                            ('_MPI',1)],
                             include_dirs=PYADH_PETSC_INCLUDE_DIRS + [numpy.get_include(),'include',
                                           PYADH_SUPERLU_INCLUDE_DIR,
					   PYADH_LAPACK_INCLUDE_DIR,
					   PYADH_BLAS_INCLUDE_DIR,
                                           PYADH_PACKAGES+"/adh/include",
                                           PYADH_PACKAGES+"/adh/main"
                                           ],
                             library_dirs=[PYADH_LAPACK_LIB_DIR,
                                           PYADH_BLAS_LIB_DIR,
                                           PYADH_PACKAGES+"/adh/lib"],
                             libraries=['m',PYADH_LAPACK_LIB,PYADH_BLAS_LIB,'adh','fileio','fric'],
                             extra_compile_args=PYADH_EXTRA_COMPILE_ARGS,
                             extra_link_args=PYADH_EXTRA_LINK_ARGS)]

    )
