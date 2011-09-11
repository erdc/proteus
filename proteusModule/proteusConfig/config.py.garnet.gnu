import os

PROTEUS_INCLUDE_DIR = os.getenv('PROTEUS_PREFIX')+'/include'
PROTEUS_LIB_DIR = os.getenv('PROTEUS_PREFIX')+'/lib'
PROTEUS_EXTRA_COMPILE_ARGS= ['-DF77_POST_UNDERSCORE','-DUSE_BLAS','-DCMRVEC_BOUNDS_CHECK','-DMV_VECTOR_BOUNDS_CHECK']
PROTEUS_EXTRA_LINK_ARGS=['-L/opt/acml/4.4.0/gnu64/lib','-lacml']
PROTEUS_EXTRA_FC_COMPILE_ARGS= ['-Wall']
PROTEUS_EXTRA_FC_LINK_ARGS=['-L/opt/acml/4.4.0/gnu64/lib','-acml']

PROTEUS_SUPERLU_INCLUDE_DIR = PROTEUS_INCLUDE_DIR
PROTEUS_SUPERLU_H   = r'"slu_ddefs.h"'
PROTEUS_SUPERLU_LIB_DIR = PROTEUS_LIB_DIR
PROTEUS_SUPERLU_LIB = 'superlu'

PROTEUS_BLAS_INCLUDE_DIR   = '.'
PROTEUS_BLAS_H     = r'"proteus_blas.h"'
PROTEUS_BLAS_LIB_DIR = '/opt/acml/4.4.0/gnu64/lib'
PROTEUS_BLAS_LIB   = 'acml'

PROTEUS_LAPACK_INCLUDE_DIR = '.'
PROTEUS_LAPACK_H   = r'"proteus_lapack.h"'
PROTEUS_LAPACK_LIB_DIR = '/opt/acml/4.4.0/gnu64/lib'
PROTEUS_LAPACK_LIB = 'acml'
PROTEUS_LAPACK_INTEGER = 'int'

PROTEUS_MPI_INCLUDE_DIR = '.'
PROTEUS_MPI_LIB_DIR = '.'
PROTEUS_MPI_LIBS = []

PROTEUS_TRIANGLE_INCLUDE_DIR = PROTEUS_INCLUDE_DIR
PROTEUS_TRIANGLE_H = r'"triangle.h"'
PROTEUS_TRIANGLE_LIB_DIR = PROTEUS_LIB_DIR
PROTEUS_TRIANGLE_LIB ='tri'


PROTEUS_DAETK_INCLUDE_DIR = [PROTEUS_INCLUDE_DIR]
PROTEUS_DAETK_LIB_DIR = PROTEUS_LIB_DIR
PROTEUS_DAETK_LIB ='daetk'
PROTEUS_PETSC_LIB_DIRS = [PROTEUS_LIB_DIR, '/opt/acml/4.4.0/gnu64/lib', '/opt/cray/mpt/5.1.2/xt/gemini/mpich2-gnu/lib', '/opt/gcc/4.5.1/snos/lib/gcc/x86_64-suse-linux/4.5.1', '/opt/gcc/4.5.1/snos/lib64', '/opt/gcc/4.5.1/snos/lib']
PROTEUS_PETSC_LIBS = ['petsc', 'X11', 'HYPRE', 'l,-rpath,/opt/gcc/4.5.1/snos/lib/gcc/x86_64-suse-linux/4.5.1', 'stdc++', 'superlu_dist_2.4', 'cmumps', 'dmumps', 'smumps', 'zmumps', 'mumps_common', 'pord', 'parmetis', 'metis', 'scalapack', 'blacs', 'superlu_4.0', 'acml', 'mpich', 'm', 'l,-rpath,/opt/gcc/4.5.1/snos/lib/gcc/x86_64-suse-linux/4.5.1', '/opt/gcc/4.5.1/snos/lib/gcc/x86_64-suse-linux/4.5.1', 'dl', 'gcc_s', 'gfortran', 'm', 'm', 'stdc++', 'dl', 'gcc_s', 'dl']
PROTEUS_PETSC_INCLUDE_DIRS = [PROTEUS_INCLUDE_DIR, '/opt/cray/mpt/5.1.2/xt/gemini/mpich2-gnu/include']
