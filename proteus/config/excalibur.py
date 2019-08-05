from __future__ import absolute_import
from .default import *

PROTEUS_PRELOAD_LIBS=[]
PROTEUS_MPI_INCLUDE_DIRS = [os.path.join(os.getenv("CRAY_MPICH2_DIR"),'include')]
PROTEUS_MPI_LIB_DIRS = [os.path.join(os.getenv("CRAY_MPICH2_DIR"),'lib')]
PROTEUS_MPI_LIBS =[]
PROTEUS_SUPERLU_LIB_DIR = os.path.join(prefix,'lib64')
