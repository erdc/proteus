#!/bin/tcsh
#
#diamond
#
module unload compilers mpi
module load mpi/sgi_mpi-2.04
module swap compiler compiler/intel11.1.072
#module load mpi/intelmpi-4.0.0
setenv PROTEUS_ARCH diamond
setenv PROTEUS_PREFIX ${PROTEUS}/${PROTEUS_ARCH}
setenv PROTEUS_PYTHON ${PROTEUS_PREFIX}/bin/python
setenv PATH .:${PROTEUS_PREFIX}/bin:${PATH}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}/lib:${PROTEUS_PREFIX}/lib:/opt/intel/cmkl/10.2.0.013/lib/em64t
