from .default import *

PROTEUS_PRELOAD_LIBS=[]
PROTEUS_EXTRA_LINK_ARGS=['-L'+os.path.join(os.getenv("CRAY_LIBSCI_PREFIX_DIR"),'lib'),'-lsci_gnu'] + platform_extra_link_args
PROTEUS_EXTRA_FC_LINK_ARGS=['-L'+os.path.join(os.getenv("CRAY_LIBSCI_PREFIX_DIR"),'lib'),'-lsci_gnu']
PROTEUS_BLAS_LIB_DIR = os.path.join(os.getenv("CRAY_LIBSCI_PREFIX_DIR"),'lib')
PROTEUS_BLAS_LIB   = 'sci_gnu'
PROTEUS_LAPACK_LIB_DIR = os.path.join(os.getenv("CRAY_LIBSCI_PREFIX_DIR"),'lib')
PROTEUS_LAPACK_LIB = 'sci_gnu'
PROTEUS_MPI_INCLUDE_DIRS = [os.path.join(os.getenv("CRAY_MPICH2_DIR"),'include')]
PROTEUS_MPI_LIB_DIRS = [os.path.join(os.getenv("CRAY_MPICH2_DIR"),'lib')]
PROTEUS_MPI_LIBS =[]
PROTEUS_SUPERLU_LIB_DIR = os.path.join(prefix,'lib64')