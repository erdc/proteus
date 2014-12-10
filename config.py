import os
from os.path import join as pjoin
import sys


prefix = os.getenv('PROTEUS_PREFIX')
if not prefix:
    prefix = sys.exec_prefix

PROTEUS_INCLUDE_DIR = pjoin(prefix, 'include')
PROTEUS_LIB_DIR = pjoin(prefix, 'lib')

platform_extra_compile_args = []
platform_extra_link_args = []
platform_blas_h = None
platform_lapack_h = None
platform_lapack_integer = None

if sys.platform == 'darwin':
    platform_extra_link_args = ['-framework', 'Accelerate']
    platform_lapack_integer = '__CLPK_integer'
    platform_blas_h = r'<veclib/cblas.h>'
    platform_lapack_h = r'<veclib/clapack.h>'
elif sys.platform == 'linux2':
    platform_extra_compile_args = ['-DPETSC_INCLUDE_AS_C']
    platform_extra_link_args = ['-Wl,-rpath,' + PROTEUS_LIB_DIR]
    platform_blas_h = r'"proteus_blas.h"'
    platform_lapack_h = r'"proteus_lapack.h"'

PROTEUS_EXTRA_COMPILE_ARGS= ['-Wall',
                             '-DF77_POST_UNDERSCORE',
                             '-DUSE_BLAS',
                             '-DCMRVEC_BOUNDS_CHECK',
                             '-DMV_VECTOR_BOUNDS_CHECK'] + platform_extra_compile_args

PROTEUS_EXTRA_LINK_ARGS=['-lblas'] + platform_extra_link_args


PROTEUS_EXTRA_FC_COMPILE_ARGS= ['-Wall']
PROTEUS_EXTRA_FC_LINK_ARGS=['-lblas']

PROTEUS_SUPERLU_INCLUDE_DIR = PROTEUS_INCLUDE_DIR
PROTEUS_SUPERLU_H   = r'"slu_ddefs.h"'
PROTEUS_SUPERLU_LIB_DIR = PROTEUS_LIB_DIR
PROTEUS_SUPERLU_LIB = 'superlu_4.1'

PROTEUS_BLAS_INCLUDE_DIR   = '.'
if platform_blas_h:
    PROTEUS_BLAS_H = platform_blas_h
else:
    PROTEUS_BLAS_H = r'"cblas.h"'
PROTEUS_BLAS_LIB_DIR = '.'
PROTEUS_BLAS_LIB   = 'blas'

PROTEUS_LAPACK_INCLUDE_DIR = '.'
if platform_lapack_h:
    PROTEUS_LAPACK_H = platform_lapack_h
else:
    PROTEUS_LAPACK_H   = r'"clapack.h"'
PROTEUS_LAPACK_LIB_DIR = '.'
PROTEUS_LAPACK_LIB = 'lapack'
if platform_lapack_integer:
    PROTEUS_LAPACK_INTEGER = platform_lapack_integer
else:
    PROTEUS_LAPACK_INTEGER = 'int'

PROTEUS_TRIANGLE_INCLUDE_DIR = PROTEUS_INCLUDE_DIR
PROTEUS_TRIANGLE_H = r'"triangle.h"'
PROTEUS_TRIANGLE_LIB_DIR = PROTEUS_LIB_DIR
PROTEUS_TRIANGLE_LIB ='tri'

PROTEUS_DAETK_INCLUDE_DIR = [PROTEUS_LIB_DIR]
PROTEUS_DAETK_LIB_DIR = PROTEUS_LIB_DIR
PROTEUS_DAETK_LIB ='daetk'

PROTEUS_MPI_INCLUDE_DIRS = []
PROTEUS_MPI_LIB_DIRS = []
PROTEUS_MPI_LIBS =[]

PROTEUS_PETSC_LIB_DIRS = []
PROTEUS_PETSC_LIBS = []
PROTEUS_PETSC_INCLUDE_DIRS = []
