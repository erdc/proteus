#
# triton intel
#
module purge
module load intel
module load openmpi_mx

export ICC_LIB=/opt/intel/Compiler/11.1/072/lib/intel64/
export MKL_LIB=/opt/intel/Compiler/11.1/072/mkl/lib/em64t/

export PROTEUS_ARCH=triton
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F77=mpif77
export F90=mpif90
export PROTEUS_PREFIX=${PROTEUS}/${PROTEUS_ARCH}
export PROTEUS_PYTHON=${PROTEUS}/${PROTEUS_ARCH}/bin/python
export PATH=${PROTEUS}/${PROTEUS_ARCH}/bin:${PATH}
export LD_LIBRARY_PATH=${PROTEUS}/${PROTEUS_ARCH}/lib:$MKL_LIB:$ICC_LIB:${LD_LIBRARY_PATH}
