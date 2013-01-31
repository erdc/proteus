#
#garnet.gnu
#
module swap PrgEnv-pgi PrgEnv-gnu
export PROTEUS_ARCH=garnet.gnu
export PROTEUS_PREFIX=${PROTEUS}/${PROTEUS_ARCH}
export PROTEUS_PYTHON=${PROTEUS_PREFIX}/bin/python
export PATH=.:${PROTEUS_PREFIX}/bin:${HOME}/src/proteus/garnet/bin:${PATH}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PROTEUS_PREFIX}/lib
export CC=gcc
export CXX=g++
export FC=gfortran
export F77=gfortran
export F90=gfortran
