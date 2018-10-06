from __future__ import absolute_import
from .default import *

PROTEUS_PRELOAD_LIBS=[]
PROTEUS_EXTRA_LINK_ARGS=['-L'+repr(os.getenv("CRAY_LIBSCI_PREFIX_DIR"))+'/lib','-lsci_gnu'] + platform_extra_link_args
PROTEUS_EXTRA_FC_LINK_ARGS=['-L'+repr(os.getenv("CRAY_LIBSCI_PREFIX_DIR"))+"/lib",'-lsci_gnu']
PROTEUS_BLAS_LIB_DIR = ''+repr(os.getenv("CRAY_LIBSCI_PREFIX_DIR"))
PROTEUS_BLAS_LIB   = 'sci_gnu'
PROTEUS_LAPACK_LIB_DIR = ''+repr(os.getenv("CRAY_LIBSCI_PREFIX_DIR"))
PROTEUS_LAPACK_LIB = 'sci_gnu'
PROTEUS_MPI_INCLUDE_DIRS = [repr(os.getenv("CRAY_MPICH2_DIR"))+'/include']
PROTEUS_MPI_LIB_DIRS = [repr(os.getenv("CRAY_MPICH2_DIR"))+'/lib']
PROTEUS_MPI_LIBS =[]
PROTEUS_SUPERLU_LIB_DIR = os.path.join(prefix,'lib64')
