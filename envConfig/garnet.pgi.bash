#
#garent.pgi
#
module swap PrgEnv-gnu PrgEnv-pgi
module load acml
module swap xt-mpich2 xt-mpich2/5.3.1
module swap xt-libsci xt-libsci/11.0.00
module swap xt-asyncpe xt-asyncpe/5.00
export PROTEUS_ARCH=garnet.pgi
export PROTEUS_PREFIX=${PROTEUS}/${PROTEUS_ARCH}
export PROTEUS_PYTHON=${PROTEUS_PREFIX}/bin/python
export PATH=.:${PROTEUS_PREFIX}/bin:${HOME}/src/proteus/garnet/bin:${PATH}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PROTEUS_PREFIX}/lib:${ACML_DIR}/pgi64/lib
