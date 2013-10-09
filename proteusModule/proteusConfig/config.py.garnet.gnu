import os

PROTEUS_INCLUDE_DIR = os.getenv('PROTEUS_PREFIX')+'/include'
PROTEUS_LIB_DIR = os.getenv('PROTEUS_PREFIX')+'/lib'
PROTEUS_EXTRA_COMPILE_ARGS= ['-DF77_POST_UNDERSCORE','-DUSE_BLAS','-DCMRVEC_BOUNDS_CHECK','-DMV_VECTOR_BOUNDS_CHECK']
PROTEUS_EXTRA_LINK_ARGS=[ '-L/opt/cray/libsci/12.1.00/gnu/48/interlagos/lib','-lsci_gnu','-lstdc++']
PROTEUS_EXTRA_FC_COMPILE_ARGS= ['-Wall']
PROTEUS_EXTRA_FC_LINK_ARGS=['-L/opt/cray/libsci/12.1.00/gnu/48/interlagos/lib','-lsci_gnu']

PROTEUS_SUPERLU_INCLUDE_DIR = PROTEUS_INCLUDE_DIR
PROTEUS_SUPERLU_H   = r'"slu_ddefs.h"'
PROTEUS_SUPERLU_LIB_DIR = PROTEUS_LIB_DIR
PROTEUS_SUPERLU_LIB = 'superlu'

PROTEUS_BLAS_INCLUDE_DIR   = '.'
PROTEUS_BLAS_H     = r'"proteus_blas.h"'
PROTEUS_BLAS_LIB_DIR = '/opt/cray/libsci/12.1.00/gnu/48/interlagos/lib'
PROTEUS_BLAS_LIB   = 'sci_gnu'

PROTEUS_LAPACK_INCLUDE_DIR = '.'
PROTEUS_LAPACK_H   = r'"proteus_lapack.h"'
PROTEUS_LAPACK_LIB_DIR = '/opt/cray/libsci/12.1.00/gnu/48/interlagos/lib'
PROTEUS_LAPACK_LIB = 'sci_gnu'
PROTEUS_LAPACK_INTEGER = 'int'

PROTEUS_MPI_INCLUDE_DIR = '.'
PROTEUS_MPI_LIB_DIR = '.'
PROTEUS_MPI_LIBS = ['mpich']

PROTEUS_TRIANGLE_INCLUDE_DIR = PROTEUS_INCLUDE_DIR
PROTEUS_TRIANGLE_H = r'"triangle.h"'
PROTEUS_TRIANGLE_LIB_DIR = PROTEUS_LIB_DIR
PROTEUS_TRIANGLE_LIB ='tri'


PROTEUS_DAETK_INCLUDE_DIR = [PROTEUS_INCLUDE_DIR]
PROTEUS_DAETK_LIB_DIR = PROTEUS_LIB_DIR
PROTEUS_DAETK_LIB ='daetk'
PROTEUS_PETSC_EXTRA_COMPILE_ARGS = []
PROTEUS_PETSC_EXTRA_LINKE_ARGS = []
PROTEUS_PETSC_LIB_DIRS = ['L/u/cekees/proteus/garnet.gnu/lib']
PROTEUS_PETSC_LIBS = ['petsc', 'parmetis', 'metis', 'pthread', 'sci', 'dl']
PROTEUS_PETSC_INCLUDE_DIRS = []
PROTEUS_PETSC_LIB_DIRS = ['L/u/cekees/proteus/garnet.gnu/lib']
PROTEUS_PETSC_LIBS = ['petsc', 'pthread', 'sci', 'dl']
PROTEUS_PETSC_INCLUDE_DIRS = []
PROTEUS_PETSC_LIB_DIRS = ['L/u/cekees/proteus/garnet.gnu/lib']
PROTEUS_PETSC_LIBS = ['petsc', 'parmetis', 'metis', 'pthread', 'sci', 'dl']
PROTEUS_PETSC_INCLUDE_DIRS = []
PROTEUS_PETSC_LIB_DIRS = ['L/u/cekees/proteus/garnet.gnu/lib', '/opt/cray/udreg/2.3.1-1.0400.4264.3.1.gem/lib64', '/opt/cray/ugni/2.3-1.0400.4374.4.88.gem/lib64', '/opt/cray/pmi/3.0.1-1.0000.8917.33.1.gem/lib64', '/opt/cray/dmapp/3.2.1-1.0400.4255.2.159.gem/lib64', '/opt/cray/xpmem/0.1-2.0400.31280.3.1.gem/lib64', '/opt/cray/mpt/6.0.0/gni/mpich2-gnu/48/lib', '/opt/cray/libsci/12.1.00/gnu/48/interlagos/lib', '/opt/fftw/3.3.0.3/interlagos/lib', '/opt/cray/xe-sysroot/4.0.46.securitypatch.20120828/usr/lib64', '/opt/cray/xe-sysroot/4.0.46.securitypatch.20120828/lib64', '/opt/cray/xe-sysroot/4.0.46.securitypatch.20120828/usr/lib/alps', '/usr/lib/alps', '/lustre/home1/opt/gcc/4.8.0/snos/lib/gcc/x86_64-suse-linux/4.8.0', '/lustre/home1/opt/gcc/4.8.0/snos/lib/gcc', '/lustre/home1/opt/gcc/4.8.0/snos/lib64', '/lustre/home1/opt/gcc/4.8.0/snos/lib', '/opt/cray/atp/1.6.3/lib']
PROTEUS_PETSC_LIBS = ['petsc', 'cmumps', 'dmumps', 'smumps', 'zmumps', 'mumps_common', 'pord', 'superlu_4.3', 'superlu_dist_3.3', 'HYPRE', 'mpichcxx_gnu_48', 'sci_gnu', 'parmetis', 'metis', 'pthread', 'gfortran', 'm', 'gfortran', 'm', 'gfortran', 'm', 'm', 'm', 'm', 'quadmath', 'm', 'mpichcxx_gnu_48', 'dl', 'AtpSigHCommData', 'AtpSigHandler', 'gfortran', 'scicpp_gnu', 'sci_gnu_mp', 'stdc++', 'fftw3_mpi', 'fftw3f_mpi', 'fftw3_threads', 'fftw3f_threads', 'fftw3', 'fftw3f', 'mpich_gnu_48', 'mpichf90_gnu_48', 'mpl', 'rt', 'xpmem', 'dmapp', 'ugni', 'pmi', 'alpslli', 'alpsutil', 'udreg', 'pthread', 'gomp', 'gcc_eh', 'dl']
PROTEUS_PETSC_INCLUDE_DIRS = []
PROTEUS_PETSC_LIB_DIRS = ['L/u/cekees/proteus/garnet.gnu/lib', '/opt/cray/udreg/2.3.1-1.0400.4264.3.1.gem/lib64', '/opt/cray/ugni/2.3-1.0400.4374.4.88.gem/lib64', '/opt/cray/pmi/3.0.1-1.0000.8917.33.1.gem/lib64', '/opt/cray/dmapp/3.2.1-1.0400.4255.2.159.gem/lib64', '/opt/cray/xpmem/0.1-2.0400.31280.3.1.gem/lib64', '/opt/cray/mpt/6.0.0/gni/mpich2-gnu/48/lib', '/opt/cray/libsci/12.1.00/gnu/48/interlagos/lib', '/opt/cray/xe-sysroot/4.0.46.securitypatch.20120828/usr/lib64', '/opt/cray/xe-sysroot/4.0.46.securitypatch.20120828/lib64', '/opt/cray/xe-sysroot/4.0.46.securitypatch.20120828/usr/lib/alps', '/usr/lib/alps', '/lustre/home1/opt/gcc/4.8.0/snos/lib/gcc/x86_64-suse-linux/4.8.0', '/lustre/home1/opt/gcc/4.8.0/snos/lib/gcc', '/lustre/home1/opt/gcc/4.8.0/snos/lib64', '/lustre/home1/opt/gcc/4.8.0/snos/lib', '/opt/cray/atp/1.6.3/lib']
PROTEUS_PETSC_LIBS = ['petsc', 'cmumps', 'dmumps', 'smumps', 'zmumps', 'mumps_common', 'pord', 'superlu_4.3', 'superlu_dist_3.3', 'HYPRE', 'mpichcxx_gnu_48', 'sci_gnu', 'parmetis', 'metis', 'pthread', 'gfortran', 'm', 'gfortran', 'm', 'gfortran', 'm', 'm', 'm', 'm', 'quadmath', 'm', 'mpichcxx_gnu_48', 'dl', 'AtpSigHCommData', 'AtpSigHandler', 'gfortran', 'scicpp_gnu', 'sci_gnu_mp', 'stdc++', 'mpich_gnu_48', 'mpichf90_gnu_48', 'mpl', 'rt', 'xpmem', 'dmapp', 'ugni', 'pmi', 'alpslli', 'alpsutil', 'udreg', 'pthread', 'gomp', 'gcc_eh', 'dl']
PROTEUS_PETSC_INCLUDE_DIRS = []
