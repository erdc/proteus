export PROTEUS_ARCH=sl6
export PROTEUS_PREFIX=${PROTEUS}/${PROTEUS_ARCH}
export PROTEUS_PYTHON=${PROTEUS_PREFIX}/bin/python
export PATH=${PROTEUS_PREFIX}/bin:${PATH}
export LD_LIBRARY_PATH=${PROTEUS_PREFIX}/lib
export DYLD_LIBRARY_PATH=${PROTEUS_PREFIX}/lib
module load gcc
module load openmpi
module load cmake
