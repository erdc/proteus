from .default import *

PROTEUS_PRELOAD_LIBS=[]
PROTEUS_EXTRA_LINK_ARGS=['-L/opt/cray/libsci/13.0.3/GNU/49/x86_64/lib','-lsci_gnu'] + platform_extra_link_args
PROTEUS_EXTRA_FC_LINK_ARGS=['-L/opt/cray/libsci/13.0.3/GNU/49/x86_64/lib','-lsci_gnu']
PROTEUS_BLAS_LIB_DIR = '/opt/cray/libsci/13.0.3/GNU/49/x86_64/lib'
PROTEUS_BLAS_LIB   = 'sci_gnu'
PROTEUS_LAPACK_LIB_DIR = '/opt/cray/libsci/13.0.3/GNU/49/x86_64/lib'
PROTEUS_LAPACK_LIB = 'sci_gnu'
PROTEUS_MPI_INCLUDE_DIRS = ['/opt/cray/mpt/7.1.3/gni/mpich2-gnu/49/include']
PROTEUS_MPI_LIB_DIRS = ['/opt/cray/mpt/7.1.3/gni/mpich2-gnu/49/lib']
PROTEUS_MPI_LIBS =[]
