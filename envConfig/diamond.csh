#!/bin/csh
#
#diamond
#
module unload compilers mpi
module load mpi/intelmpi-4.0.0
setenv PROTEUS_ARCH diamond
setenv PROTEUS_PREFIX ${PROTEUS}/${PROTEUS_ARCH}
setenv PROTEUS_PYTHON ${PROTEUS_PREFIX}/bin/python
setenv CC mpicc
setenv CXX mpicxx
setenv FC mpifc
setenv F77 mpif77
setenv F90 mpif90
setenv CFLAGS "-fPIC"
setenv CXXFLAGS "-fPIC"
setenv FCFLAGS "-fPIC"
setenv F77FLAGS "-fPIC"
setenv F90FLAGS "-fPIC"
setenv FFLAGS "-fPIC"
setenv LD "mpicxx"
setenv LDFLAGS "-lirc"
setenv PATH .:${PROTEUS_PREFIX}/bin:${PATH}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PROTEUS_PREFIX}/lib
