from .default import *

platform_extra_link_args = ['-Wl,'+ "-L" + PROTEUS_LIB_DIR]
PROTEUS_EXTRA_LINK_ARGS=['-lopenblas'] + platform_extra_link_args
PROTEUS_EXTRA_FC_LINK_ARGS=['-lopenblas']
PROTEUS_BLAS_INCLUDE_DIR   = prefix + "/include"
PROTEUS_BLAS_LIB   = 'openblas'
PROTEUS_LAPACK_INCLUDE_DIR = prefix + "/include"
PROTEUS_MPI_INCLUDE_DIRS = ["/opt/apps/gcc4_9/mvapich2/2.1" + "/include"]
PROTEUS_MPI_LIB_DIRS = ["/opt/apps/gcc4_9/mvapich2/2.1" + "/lib"]
PROTEUS_MPI_LIBS =[]
