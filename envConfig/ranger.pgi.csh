#!/bin/csh
#
#ranger.pgi
#
module load acml
module load git
setenv PROTEUS_ARCH ranger.pgi
setenv CC  pgcc
setenv CXX pgCC
setenv FC  pgf90
setenv F77 pgf77
setenv F90 pgf90
#
setenv PROTEUS_PREFIX ${PROTEUS}/${PROTEUS_ARCH}
setenv PROTEUS_PYTHON ${PROTEUS}/${PROTEUS_ARCH}/bin/python
setenv PATH ${PROTEUS}/${PROTEUS_ARCH}/bin:${PATH}
setenv LD_LIBRARY_PATH ${PROTEUS}/${PROTEUS_ARCH}/lib:${LD_LIBRARY_PATH}
