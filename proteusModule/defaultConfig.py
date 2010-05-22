import os

PYADH_PACKAGES = os.getenv('PYADH_PACKAGES',os.getenv('HOME')+'/src/pyadh-packages')
PYADH_EXTRA_COMPILE_ARGS= ['-Wall']
PYADH_EXTRA_LINK_ARGS=['-lblas','-lstdc++']

PYADH_SUPERLU_INCLUDE_DIR = PYADH_PACKAGES+'/SuperLU/SRC'
PYADH_SUPERLU_H   = r'"slu_ddefs.h"'
PYADH_SUPERLU_LIB_DIR = PYADH_PACKAGES+'/SuperLU'
PYADH_SUPERLU_LIB = 'superlu'

PYADH_BLAS_INCLUDE_DIR   = '.'
PYADH_BLAS_H     = r'"vecLib/cblas.h"'
PYADH_BLAS_LIB_DIR = '.'
PYADH_BLAS_LIB   = 'm'

PYADH_LAPACK_INCLUDE_DIR = '.'
PYADH_LAPACK_H   = r'"vecLib/clapack.h"'
PYADH_LAPACK_LIB_DIR = '.'
PYADH_LAPACK_LIB = 'm'
PYADH_LAPACK_INTEGER = '__CLPK_integer'

PYADH_TRIANGLE_INCLUDE_DIR = PYADH_PACKAGES+'/triangle/src'
PYADH_TRIANGLE_H = r'"triangle.h"'
PYADH_TRIANGLE_LIB_DIR = PYADH_PACKAGES+'/triangle/lib'
PYADH_TRIANGLE_LIB ='tri'

PYADH_DAETK_INCLUDE_DIR = [PYADH_PACKAGES+'/daetk',PYADH_PACKAGES+'/daetk/pete/pete-2.1.0/src']
PYADH_DAETK_LIB_DIR = PYADH_PACKAGES+'/daetk'
PYADH_DAETK_LIB ='daetk'

PYADH_PETSC_INCLUDE_DIRS = [PYADH_PACKAGES+'/petsc/include',
                            PYADH_PACKAGES+'/petsc/bmake/'+os.getenv('PETSC_ARCH')]
PYADH_PETSC_LIB_DIRS = [PYADH_PACKAGES+'/petsc/lib/'+os.getenv('PETSC_ARCH'),
                        PYADH_PACKAGES+'/petsc/externalpackages/ParMetis-dev-p1/'+os.getenv('PETSC_ARCH')+'/lib',
                        '/usr/X11R6/lib']
PYADH_PETSC_LIBS =['petsccontrib','petscts','petscsnes','petscksp','petscdm','petscmat','petscvec','petsc','parmetis','metis','X11']

PYADH_MPI_INCLUDE_DIR = '.'
PYADH_MPI_LIB_DIR = '.'
PYADH_MPI_LIBS = []
#PYADH_MPI_INCLUDE_DIR = '/sw/include'
#PYADH_MPI_LIB_DIR = '/sw/lib/openmpi'
#PYADH_MPI_LIBS =['mpi','mpi_f77','mpi_f90','mpi_cxx','open-rte','open-pal']
# PYADH_MPI_INCLUDE_DIR = '/usr/local/include'
# PYADH_MPI_LIB_DIR = '/usr/local/lib'
# PYADH_MPI_LIBS =['mpi','orte','opal','dl']

