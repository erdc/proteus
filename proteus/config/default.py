import os
from os.path import join as pjoin
import sys
import platform
PROTEUS_PRELOAD_LIBS=[]

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
    platform_blas_h = r'<Accelerate/Accelerate.h>'
    platform_lapack_h = r'<Accelerate/Accelerate.h>'
    platform_extra_compile_args = ['-flax-vector-conversions','-DPETSC_INCLUDE_AS_C','-DPETSC_SKIP_COMPLEX']
    major,minor = platform.mac_ver()[0].split('.')[0:2]
    os.environ["MACOSX_DEPLOYMENT_TARGET"]= major+'.'+minor
elif sys.platform == 'linux2':
    platform_extra_compile_args = ['-DPETSC_INCLUDE_AS_C', '-DPETSC_SKIP_COMPLEX']
    platform_extra_link_args = ['-L'+PROTEUS_LIB_DIR,'-Wl,-rpath,' + PROTEUS_LIB_DIR]
    platform_blas_h = r'"proteus_blas.h"'
    platform_lapack_h = r'"proteus_lapack.h"'

PROTEUS_EXTRA_COMPILE_ARGS= ['-Wall',
                             '-DF77_POST_UNDERSCORE',
                             '-DUSE_BLAS',
                             '-DCMRVEC_BOUNDS_CHECK',
                             '-DMV_VECTOR_BOUNDS_CHECK'] + platform_extra_compile_args

def get_flags(package):
    """ Checks the environment for presence of PACKAGE_DIR

    And either returns PACKAGE_DIR/[include, lib] or the Proteus include flags.

    This supports building Proteus using packages provides via environment variables.
    """

    package_dir_env = os.getenv(package.upper() + '_DIR')
    if package_dir_env:
        include_dir = pjoin(package_dir_env, 'include')
        lib_dir = pjoin(package_dir_env, 'lib')
    else:
        include_dir = PROTEUS_INCLUDE_DIR
        lib_dir = PROTEUS_LIB_DIR
    return include_dir, lib_dir

PROTEUS_BLAS_INCLUDE_DIR, PROTEUS_BLAS_LIB_DIR = get_flags('blas')

if sys.platform == 'darwin':
    PROTEUS_BLAS_LIB ='m'
    PROTEUS_BLAS_LIB_DIR = PROTEUS_LIB_DIR
    PROTEUS_BLAS_INCLUDE_DIR = PROTEUS_INCLUDE_DIR
    PROTEUS_EXTRA_LINK_ARGS=platform_extra_link_args
else:
    PROTEUS_BLAS_LIB   ='openblas'
    PROTEUS_BLAS_INCLUDE_DIR, PROTEUS_BLAS_LIB_DIR = get_flags('blas')

PROTEUS_EXTRA_LINK_ARGS=[]

PROTEUS_CHRONO_INCLUDE_DIR, PROTEUS_CHRONO_LIB_DIR = get_flags('chrono')

PROTEUS_EXTRA_FC_COMPILE_ARGS= ['-Wall']
PROTEUS_EXTRA_FC_LINK_ARGS=[]



if platform_blas_h:
    PROTEUS_BLAS_H = platform_blas_h
else:
    PROTEUS_BLAS_H = r'"cblas.h"'

PROTEUS_LAPACK_INCLUDE_DIR, PROTEUS_LAPACK_LIB_DIR = get_flags('lapack')
if platform_lapack_h:
    PROTEUS_LAPACK_H = platform_lapack_h
else:
    PROTEUS_LAPACK_H   = r'"clapack.h"'

if sys.platform == 'darwin':
    PROTEUS_LAPACK_LIB ='m'
else:
    PROTEUS_LAPACK_LIB = 'openblas'

if platform_lapack_integer:
    PROTEUS_LAPACK_INTEGER = platform_lapack_integer
else:
    PROTEUS_LAPACK_INTEGER = 'int'


PROTEUS_TRIANGLE_INCLUDE_DIR, PROTEUS_TRIANGLE_LIB_DIR = get_flags('triangle')
PROTEUS_TRIANGLE_H = r'"triangle.h"'
PROTEUS_TRIANGLE_LIB ='tri'

PROTEUS_DAETK_INCLUDE_DIR, PROTEUS_DAETK_LIB_DIR = get_flags('daetk')
PROTEUS_DAETK_LIB ='daetk'
PROTEUS_DAETK_LIB_DIRS = [PROTEUS_DAETK_LIB_DIR]

PROTEUS_MPI_INCLUDE_DIR, PROTEUS_MPI_LIB_DIR = get_flags('mpi')
PROTEUS_MPI_INCLUDE_DIRS = [PROTEUS_MPI_INCLUDE_DIR]
PROTEUS_MPI_LIB_DIRS = [PROTEUS_MPI_LIB_DIR]
PROTEUS_MPI_LIBS =[]

PROTEUS_PETSC_INCLUDE_DIR, PROTEUS_PETSC_LIB_DIR = get_flags('petsc')
PROTEUS_PETSC_LIB_DIRS = [PROTEUS_PETSC_LIB_DIR]
PROTEUS_PETSC_LIBS = ['petsc']
PROTEUS_PETSC_INCLUDE_DIRS = [PROTEUS_PETSC_INCLUDE_DIR]

PROTEUS_SCOREC_INCLUDE_DIR, PROTEUS_SCOREC_LIB_DIR = get_flags('scorec')
PROTEUS_PARMETIS_INCLUDE_DIR, PROTEUS_PARMETIS_LIB_DIR = get_flags('parmetis')
PROTEUS_ZOLTAN_INCLUDE_DIR, PROTEUS_ZOLTAN_LIB_DIR = get_flags('zoltan')
PROTEUS_SCOREC_INCLUDE_DIRS = [PROTEUS_SCOREC_INCLUDE_DIR, PROTEUS_PETSC_INCLUDE_DIR, PROTEUS_ZOLTAN_INCLUDE_DIR, PROTEUS_PARMETIS_INCLUDE_DIR, PROTEUS_MPI_INCLUDE_DIR, PROTEUS_LAPACK_INCLUDE_DIR, PROTEUS_BLAS_INCLUDE_DIR]
PROTEUS_SCOREC_LIB_DIRS =     [PROTEUS_SCOREC_LIB_DIR,     PROTEUS_PETSC_LIB_DIR,     PROTEUS_ZOLTAN_LIB_DIR, PROTEUS_PARMETIS_LIB_DIR,     PROTEUS_MPI_LIB_DIR, PROTEUS_LAPACK_LIB_DIR, PROTEUS_BLAS_LIB_DIR]
PROTEUS_SCOREC_LIBS = [
    'spr',
    'ma',
    'parma',
    'apf_zoltan',
    'mds',
    'apf',
    'mth',
    'gmi',
    'pcu',
    'lion',
    'zoltan',
    'parmetis',
    'metis',
    'sam']+PROTEUS_PETSC_LIBS+PROTEUS_MPI_LIBS

PROTEUS_SCOREC_EXTRA_LINK_ARGS = []
PROTEUS_SCOREC_EXTRA_COMPILE_ARGS = ['-g','-DMESH_INFO']

if os.getenv('SIM_INCLUDE_DIR') is not None:
  PROTEUS_SCOREC_INCLUDE_DIRS = PROTEUS_SCOREC_INCLUDE_DIRS+[pjoin(prefix, 'include'), os.getenv('SIM_INCLUDE_DIR')]
  PROTEUS_SCOREC_LIB_DIRS = PROTEUS_SCOREC_LIB_DIRS+[pjoin(prefix, 'lib'), os.getenv('SIM_LIB_DIR')]
  PROTEUS_SCOREC_LIBS = PROTEUS_SCOREC_LIBS
  PROTEUS_SCOREC_EXTRA_COMPILE_ARGS = PROTEUS_SCOREC_EXTRA_COMPILE_ARGS

PROTEUS_SUPERLU_INCLUDE_DIR = PROTEUS_PETSC_INCLUDE_DIR
PROTEUS_SUPERLU_LIB_DIR = PROTEUS_PETSC_LIB_DIR
PROTEUS_SUPERLU_H   = r'"slu_ddefs.h"'
PROTEUS_SUPERLU_LIB = 'superlu'

PROTEUS_HDF5_INCLUDE_DIR, PROTEUS_HDF5_LIB_DIR = get_flags('hdf5')
PROTEUS_HDF5_LIB_DIRS = [PROTEUS_HDF5_LIB_DIR]
PROTEUS_HDF5_LIBS = ['hdf5']
PROTEUS_HDF5_INCLUDE_DIRS = [PROTEUS_HDF5_INCLUDE_DIR]

if sys.platform == 'linux2':
    PROTEUS_EXTRA_LINK_ARGS += ['-Wl,-rpath,' + PROTEUS_HDF5_LIB_DIR]
    PROTEUS_EXTRA_LINK_ARGS += ['-Wl,-rpath,' + PROTEUS_PETSC_LIB_DIR]
    PROTEUS_EXTRA_LINK_ARGS += ['-Wl,-rpath,' + PROTEUS_SCOREC_LIB_DIR]
    PROTEUS_EXTRA_LINK_ARGS += ['-Wl,-rpath,' + PROTEUS_PARMETIS_LIB_DIR]
    PROTEUS_EXTRA_LINK_ARGS += ['-Wl,-rpath,' + PROTEUS_ZOLTAN_LIB_DIR]
    PROTEUS_EXTRA_LINK_ARGS += ['-Wl,-rpath,' + PROTEUS_LAPACK_LIB_DIR]
    PROTEUS_EXTRA_LINK_ARGS += ['-Wl,-rpath,' + PROTEUS_BLAS_LIB_DIR]
