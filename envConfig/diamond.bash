#!/bin/tcsh
#
#diamond
#
module unload compiler mpi
module load compiler/intel12.1.003
module load mpi/intelmpi-4.0.3
export PROTEUS_ARCH=diamond
export PROTEUS_PREFIX=${PROTEUS}/${PROTEUS_ARCH}
export PROTEUS_PYTHON=${PROTEUS_PREFIX}/bin/python
export PATH=${PROTEUS_PREFIX}/bin:${PATH}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}/lib:${PROTEUS_PREFIX}/lib:/opt/intel/cmkl/10.2.4.032/lib/em64t:/usr/local/applic/intel/mkl/10.2.4.032/lib/em64t:/opt/intel/Compiler/12.1.003/compiler/lib/intel64/
unset CC
unset CXX
unset F77
unset F90
unset FC
unset LD
