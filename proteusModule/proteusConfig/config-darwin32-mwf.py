import os

PYADH_PACKAGES = os.getenv('PYADH_PACKAGES',os.getenv('HOME')+'/src/proteus-packages')
PYADH_EXTRA_COMPILE_ARGS= ['-D_REENTRANT']
PYADH_EXTRA_LINK_ARGS=[]
PYADH_SUPERLU_INCLUDE_DIR = PYADH_PACKAGES+'/SuperLU/SRC'
PYADH_SUPERLU_H   = '\"slu_ddefs.h\"'
PYADH_SUPERLU_LIB_DIR = PYADH_PACKAGES+'/SuperLU'
PYADH_SUPERLU_LIB = 'superlu'

PYADH_BLAS_INCLUDE_DIR   = '.'
PYADH_BLAS_H     = '\"vecLib/cblas.h\"'
PYADH_BLAS_LIB_DIR = '.'
PYADH_BLAS_LIB   = 'm'

PYADH_LAPACK_INCLUDE_DIR = '.'
PYADH_LAPACK_H   = '\"vecLib/clapack.h\"'
PYADH_LAPACK_LIB_DIR = '.'
PYADH_LAPACK_LIB = 'm'
PYADH_LAPACK_INTEGER = '__CLPK_integer'

PYADH_TRIANGLE_INCLUDE_DIR = PYADH_PACKAGES+'/triangle/src'
PYADH_TRIANGLE_H = '\"triangle.h\"'
PYADH_TRIANGLE_LIB_DIR = PYADH_PACKAGES+'/triangle/lib'
PYADH_TRIANGLE_LIB ='tri'

#mwf added
PYADH_PETSC_DIR = os.getenv('PETSC_DIR',os.getenv('HOME')+'/src/petsc')
PYADH_PETSC_ARCH = os.getenv('PETSC_ARCH','darwin8.11.1-c-debug')

PYADH_DAETK_ARCH = os.getenv('PETSC_ARCH','darwin')
PYADH_DAETK_DIR  = os.getenv('DAETK_DIR',os.getenv('HOME')+'/src/daetk')

PYADH_DAETK_INCLUDE_DIR = [PYADH_DAETK_DIR,PYADH_DAETK_DIR+'/pete/pete-2.1.0/src']
PYADH_DAETK_LIB_DIR = PYADH_DAETK_DIR
PYADH_DAETK_LIB ='daetk'

PYADH_PETSC_INCLUDE_DIR = PYADH_PETSC_DIR+'/include'
PYADH_PETSC_LIB_DIRS = [PYADH_PETSC_DIR+'/lib/'+PYADH_PETSC_ARCH,
                        '/usr/X11R6/lib']
#PYADH_PETSC_DIR+'/externalpackages/ParMetis-dev-p1',
PYADH_PETSC_LIBS =['petscts','petscsnes','petscksp','petscdm','petscmat','petscvec','petsc','X11']
#,'metis','parmetis'

PYADH_MPI_INCLUDE_DIR = '/sw/include'
PYADH_MPI_LIB_DIR = ['/sw/lib','/sw/lib/openmpi','/sw/lib/gcc4.2/lib']
PYADH_MPI_LIBS =['mpi','orte','opal','dl','mpi_f77','mpi_f90']





