#!/bin/tcsh
#
#diamond
#
module unload compilers mpi
module load compiler/intel11.1.072
module load mpi/intelmpi-4.0.0
setenv PROTEUS_ARCH diamond
setenv PROTEUS_PREFIX ${PROTEUS}/${PROTEUS_ARCH}
setenv PROTEUS_PYTHON ${PROTEUS_PREFIX}/bin/python
setenv PATH .:${PROTEUS_PREFIX}/bin:${PATH}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}/lib:${PROTEUS_PREFIX}/lib:/usr/local/applic/intel/mkl/10.2.4.032
