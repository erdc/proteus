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
    ext_package='proteus',
    ext_modules = [Extension("waveFunctions",['proteus/waveFunctions.pyx','proteus/transportCoefficients.c'],
                            include_dirs=[numpy.get_include(),'include']),
                   Extension("subsurfaceTransportFunctions",['proteus/subsurfaceTransportFunctions.pyx'],
                            include_dirs=[numpy.get_include(),'include']),
                   Extension("pskRelations",['proteus/pskRelations.pyx'],
                             include_dirs=[numpy.get_include(),'include']),
                   Extension('cadh',
                             ['proteus/cadh.pyx','proteus/cadhimpl.c'],
                             define_macros=[('PROTEUS_SUPERLU_H',PROTEUS_SUPERLU_H),
                                            ('PROTEUS_LAPACK_H',PROTEUS_LAPACK_H),
                                            ('PROTEUS_LAPACK_INTEGER',PROTEUS_LAPACK_INTEGER),
                                            ('PROTEUS_BLAS_H',PROTEUS_BLAS_H),
                                            ('_DEBUG',1),
                                            ('_MESSG',1),
                                            ('_MPI',1)],
                             include_dirs=PROTEUS_PETSC_INCLUDE_DIRS + [numpy.get_include(),'include',
                                           PROTEUS_SUPERLU_INCLUDE_DIR,
					   PROTEUS_LAPACK_INCLUDE_DIR,
					   PROTEUS_BLAS_INCLUDE_DIR,
                                           PROTEUS_PACKAGES+"/adh/include",
                                           PROTEUS_PACKAGES+"/adh/main"
                                           ],
                             library_dirs=[PROTEUS_LAPACK_LIB_DIR,
                                           PROTEUS_BLAS_LIB_DIR,
                                           PROTEUS_PACKAGES+"/adh/lib"],
                             libraries=['m',PROTEUS_LAPACK_LIB,PROTEUS_BLAS_LIB,'adh','fileio','fric'],
                             extra_compile_args=PROTEUS_EXTRA_COMPILE_ARGS,
                             extra_link_args=PROTEUS_EXTRA_LINK_ARGS)]

    )
