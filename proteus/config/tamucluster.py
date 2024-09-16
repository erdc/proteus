from .default import *

PROTEUS_MPI_INCLUDE_DIR, PROTEUS_MPI_LIB_DIR = get_flags('mpi')
PROTEUS_MPI_INCLUDE_DIRS = [PROTEUS_MPI_INCLUDE_DIR,'/apps/openmpi/1.6.5/include']
PROTEUS_MPI_LIB_DIRS = [PROTEUS_MPI_LIB_DIR,'/apps/openmpi/1.6.5/lib64']
PROTEUS_MPI_LIBS =[]
