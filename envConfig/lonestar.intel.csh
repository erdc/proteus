#
#intel
#
module load mkl
setenv PROTEUS_ARCH lonestar.intel
setenv PROTEUS_PREFIX ${PROTEUS}/${PROTEUS_ARCH}
setenv PROTEUS_PYTHON ${PROTEUS_PREFIX}/bin/python
setenv PATH .:${PROTEUS_PREFIX}/bin:${HOME}/bin:${PATH}
setenv MV2_ON_DEMAND_THRESHOLD 2048
