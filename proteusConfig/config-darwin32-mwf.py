import os

PROTEUS_PACKAGES = os.getenv('PROTEUS_PACKAGES',os.getenv('HOME')+'/src/proteus-packages')
PROTEUS_EXTRA_COMPILE_ARGS= ['-D_REENTRANT']
PROTEUS_EXTRA_LINK_ARGS=[]
PROTEUS_SUPERLU_INCLUDE_DIR = PROTEUS_PACKAGES+'/SuperLU/SRC'
PROTEUS_SUPERLU_H   = '\"slu_ddefs.h\"'
PROTEUS_SUPERLU_LIB_DIR = PROTEUS_PACKAGES+'/SuperLU'
PROTEUS_SUPERLU_LIB = 'superlu'

PROTEUS_BLAS_INCLUDE_DIR   = '.'
PROTEUS_BLAS_H     = '\"vecLib/cblas.h\"'
PROTEUS_BLAS_LIB_DIR = '.'
PROTEUS_BLAS_LIB   = 'm'

PROTEUS_LAPACK_INCLUDE_DIR = '.'
PROTEUS_LAPACK_H   = '\"vecLib/clapack.h\"'
PROTEUS_LAPACK_LIB_DIR = '.'
PROTEUS_LAPACK_LIB = 'm'
PROTEUS_LAPACK_INTEGER = '__CLPK_integer'

PROTEUS_TRIANGLE_INCLUDE_DIR = PROTEUS_PACKAGES+'/triangle/src'
PROTEUS_TRIANGLE_H = '\"triangle.h\"'
PROTEUS_TRIANGLE_LIB_DIR = PROTEUS_PACKAGES+'/triangle/lib'
PROTEUS_TRIANGLE_LIB ='tri'

#mwf added
PROTEUS_PETSC_DIR = os.getenv('PETSC_DIR',os.getenv('HOME')+'/src/petsc')
PROTEUS_PETSC_ARCH = os.getenv('PETSC_ARCH','darwin8.11.1-c-debug')

PROTEUS_DAETK_ARCH = os.getenv('PETSC_ARCH','darwin')
PROTEUS_DAETK_DIR  = os.getenv('DAETK_DIR',os.getenv('HOME')+'/src/daetk')

PROTEUS_DAETK_INCLUDE_DIR = [PROTEUS_DAETK_DIR,PROTEUS_DAETK_DIR+'/pete/pete-2.1.0/src']
PROTEUS_DAETK_LIB_DIR = PROTEUS_DAETK_DIR
PROTEUS_DAETK_LIB ='daetk'

PROTEUS_PETSC_INCLUDE_DIR = PROTEUS_PETSC_DIR+'/include'
PROTEUS_PETSC_LIB_DIRS = [PROTEUS_PETSC_DIR+'/lib/'+PROTEUS_PETSC_ARCH,
                        '/usr/X11R6/lib']
#PROTEUS_PETSC_DIR+'/externalpackages/ParMetis-dev-p1',
PROTEUS_PETSC_LIBS =['petscts','petscsnes','petscksp','petscdm','petscmat','petscvec','petsc','X11']
#,'metis','parmetis'

PROTEUS_MPI_INCLUDE_DIR = '/sw/include'
PROTEUS_MPI_LIB_DIR = ['/sw/lib','/sw/lib/openmpi','/sw/lib/gcc4.2/lib']
PROTEUS_MPI_LIBS =['mpi','orte','opal','dl','mpi_f77','mpi_f90']
