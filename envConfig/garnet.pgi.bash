#
#garent.pgi
#
module load acml
export PROTEUS_ARCH=garnet.pgi
export PROTEUS_PREFIX=${PROTEUS}/${PROTEUS_ARCH}
export PROTEUS_PYTHON=${PROTEUS_PREFIX}/bin/python
export PATH=.:${PROTEUS_PREFIX}/bin:${HOME}/src/proteus/garnet/bin:${PATH}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PROTEUS_PREFIX}/lib:${ACML_DIR}/pgi64/lib
export CC=pgcc
export CXX=pgCC
export FC=pgf90
export F77=pgf77
export F90=pgf90
