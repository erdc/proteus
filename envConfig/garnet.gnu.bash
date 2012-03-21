#
#garnet.gnu
#
module swap PrgEnv-pgi PrgEnv-gnu
#module load acml
#module swap xt-mpich2 xt-mpich2/5.3.1
#module swap gcc gcc/4.6.0
#module swap xt-libsci xt-libsci/11.0.00
#module swap xt-asyncpe xt-asyncpe/5.00
export PROTEUS_ARCH=garnet.gnu
export PROTEUS_PREFIX=${PROTEUS}/${PROTEUS_ARCH}
export PROTEUS_PYTHON=${PROTEUS_PREFIX}/bin/python
export PATH=.:${PROTEUS_PREFIX}/bin:${HOME}/src/proteus/garnet/bin:${PATH}
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PROTEUS_PREFIX}/lib:${ACML_DIR}/gnu64/lib
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PROTEUS_PREFIX}/lib:${ACML_DIR}/gnu64/lib:/opt/xt-libsci/11.0.0/gnu/46/mc8/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PROTEUS_PREFIX}/lib:/opt/xt-libsci/10.4.9/gnu/lib/45
