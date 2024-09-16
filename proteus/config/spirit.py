from .default import *

PROTEUS_PRELOAD_LIBS=[]

PROTEUS_EXTRA_LINK_ARGS=[
'-L/app/intel/compilers/15/mkl/lib/intel64',
'-lmkl_rt',
'-L/app/intel/compilers/15/lib/intel64','-lsvml','-lirc','-limf',
'-lpthread','-lm',
'-L/app/gmpapp/gcc/platform/gcc-4.8.4/lib64',
'-lgfortran',
] + platform_extra_link_args
PROTEUS_EXTRA_FC_LINK_ARGS=['-L/app/intel/comilers/15/mkl/lib/intel64','-lmkl_rt']
PROTEUS_BLAS_LIB_DIR = '/app/intel/compilers/15/mkl/lib/intel64/'
PROTEUS_BLAS_LIB   = 'mkl_rt'
PROTEUS_LAPACK_LIB_DIR = '/app/intel/compilers/15/mkl/lib/intel64/'
PROTEUS_LAPACK_LIB = 'mkl_rt'
