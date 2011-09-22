#
#ranger.intel
#
module swap pgi intel
#module swap mvapich openmpi
module load mkl
setenv PROTEUS_ARCH ranger.intel
setenv CC mpicc
setenv CXX mpicxx
setenv FC mpif90
setenv F77 mpif77
setenv F90 mpif90
setenv PROTEUS_PREFIX ${PROTEUS}/${PROTEUS_ARCH}
setenv PROTEUS_PYTHON ${PROTEUS}/${PROTEUS_ARCH}/bin/python
setenv PATH ${PROTEUS}/${PROTEUS_ARCH}/bin:${PATH}
setenv LD_LIBRARY_PATH ${PROTEUS}/${PROTEUS_ARCH}/lib:${TACC_MKL_LIB}:${LD_LIBRARY_PATH}
