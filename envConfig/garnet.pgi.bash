#
#garent.pgi
#
module swap PrgEnv-gnu PrgEnv-pgi
module load gcc
export PROTEUS_ARCH=garnet.pgi
export PROTEUS_PREFIX=${PROTEUS}/${PROTEUS_ARCH}
export PROTEUS_PYTHON=${PROTEUS_PREFIX}/bin/python
export PATH=.:${PROTEUS_PREFIX}/bin:${HOME}/src/proteus/garnet/bin:${PATH}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PROTEUS_PREFIX}/lib:/opt/xt-libsci/10.4.9/pgi/lib
