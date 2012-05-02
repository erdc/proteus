source /app/modules/init/csh

module swap compiler/pgi/10.9  compiler/gcc/4.1
module swap mpi/pgi/openmpi/1.4.3 mpi/gnu/openmpi/1.4.3
#module load acml

setenv PROTEUS ${HOME}/src/proteus
setenv PROTEUS_ARCH hpc_util
setenv PROTEUS_PREFIX ${PROTEUS}/${PROTEUS_ARCH}
setenv PROTEUS_PYTHON ${PROTEUS_PREFIX}/bin/python
setenv PATH ${PROTEUS_PREFIX}/bin:${PATH}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PROTEUS_PREFIX}/lib
