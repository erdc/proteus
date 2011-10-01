# trestles at sdsc
#
# NOTE: There is an error in the auto-generated proteusModule/config.py
# QaD FIX: Manually change lib ==> include  in the last item of the last line
#

module purge
module load intel
module load openmpi

export ICC_LIB=/opt/intel/Compiler/11.1/072/lib/intel64/
export MKL_LIB=/opt/intel/Compiler/11.1/072/mkl/lib/em64t/

export PROTEUS_ARCH=trestles
export PROTEUS_PREFIX=${PROTEUS}/${PROTEUS_ARCH}
export PROTEUS_PYTHON=${PROTEUS_PREFIX}/bin/python
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PROTEUS_PREFIX}/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:${PROTEUS_PREFIX}/lib
export PATH=${PROTEUS_PREFIX}/bin:${PATH}
