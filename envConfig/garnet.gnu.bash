#
#garnet.gnu
#
module swap PrgEnv-pgi PrgEnv-gnu
export PROTEUS_ARCH=garnet.gnu
export PROTEUS_PREFIX=${PROTEUS}/${PROTEUS_ARCH}
export PROTEUS_PYTHON=${PROTEUS_PREFIX}/bin/python
export PATH=.:${PROTEUS_PREFIX}/bin:${HOME}/src/proteus/garnet/bin:${PATH}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PROTEUS_PREFIX}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PROTEUS_PREFIX}/lib:/opt/cray/libsci/12.1.00/gnu/48/interlagos/lib:/opt/cray/mpt/6.0.0/gni/mpich2-gnu/48/lib:/opt/fftw/3.3.0.3/interlagos/lib:/usr/local/usp/openssl/lib
export CC=gcc
export CXX=g++
export FC=gfortran
export F77=gfortran
export F90=gfortran
