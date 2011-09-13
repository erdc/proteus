#
#ranger.intel
#
module swap pgi intel
module swap mvapich openmpi
module load atlas
setenv PROTEUS_ARCH ranger.intel
setenv PATH ${TACC_ATLAS_LIB}:${PATH}
#
setenv PROTEUS_PREFIX ${PROTEUS}/${PROTEUS_ARCH}
setenv PROTEUS_PYTHON ${PROTEUS}/${PROTEUS_ARCH}/bin/python
setenv PATH ${PROTEUS}/${PROTEUS_ARCH}/bin:${PATH}
setenv LD_LIBRARY_PATH ${PROTEUS}/${PROTEUS_ARCH}/lib:${LD_LIBRARY_PATH}
