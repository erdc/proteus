from default import *

PROTEUS_PRELOAD_LIBS=[]
PROTEUS_EXTRA_LINK_ARGS=[
'-L/app/intel/compilers/14/mkl/lib/intel64',
'-lmkl_rt',
'-L/app/intel/compilers/14/lib/intel64','-lsvml','-lirc','-limf',
'-lpthread','-lm'
] + platform_extra_link_args
PROTEUS_EXTRA_FC_LINK_ARGS=['-L/app/intel/comilers/14/mkl/lib/intel64','-lmkl_rt']
PROTEUS_BLAS_LIB_DIR = '/app/intel/compilers/14/mkl/lib/intel64/'
PROTEUS_BLAS_LIB   = 'mkl_rt'
PROTEUS_LAPACK_LIB_DIR = '/app/intel/compilers/14/mkl/lib/intel64/'
PROTEUS_LAPACK_LIB = 'mkl_rt'
