#
#intel
#
module load mkl
#module swap mvapich2 openmpi
setenv CC  icc
setenv CFLAGS '-fPIC'
setenv CXX icpc
setenv CXXFLAGS '-fPIC'
setenv FC  ifort
setenv FCFLAGS '-fPIC'
setenv FFLAGS '-fPIC'
setenv F77 ifort
setenv F77FLAGS '-fPIC'
setenv F90 ifort
setenv F90FLAGS '-fPIC'
setenv PROTEUS_ARCH lonestar.intel
setenv PROTEUS_PREFIX ${PROTEUS}/${PROTEUS_ARCH}
setenv PROTEUS_PYTHON ${PROTEUS_PREFIX}/bin/python
setenv PATH .:${PROTEUS_PREFIX}/bin:${HOME}/bin:${PATH}
setenv LD_LIBRARY_PATH ${PROTEUS_PREFIX}/lib:${TACC_MKL_LIB}:/opt/apps/limic2/0.5.4/lib:${LD_LIBRARY_PATH}
