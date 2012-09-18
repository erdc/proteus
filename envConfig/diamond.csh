#!/bin/tcsh
#
#diamond
#
module unload compiler mpi
module load compiler/intel12.1.003
module load mpi/intelmpi-4.0.3
setenv PROTEUS_ARCH diamond
setenv PROTEUS_PREFIX ${PROTEUS}/${PROTEUS_ARCH}
setenv PROTEUS_PYTHON ${PROTEUS_PREFIX}/bin/python
setenv PATH .:${PROTEUS_PREFIX}/bin:${PATH}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}/lib:${PROTEUS_PREFIX}/lib:/opt/intel/cmkl/10.2.4.032/lib/em64t:/usr/local/applic/intel/mkl/10.2.4.032/lib/em64t:/opt/intel/Compiler/12.1.003/compiler/lib/intel64/
unsetenv CC
unsetenv CXX
unsetenv F77
unsetenv F90
unsetenv FC
unsetenv LD
