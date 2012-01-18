#
#gnu
#
module swap intel gcc
module load atlas
setenv PROTEUS_ARCH lonestar.gnu
setenv PROTEUS_PREFIX ${PROTEUS}/${PROTEUS_ARCH}
setenv PROTEUS_PYTHON ${PROTEUS_PREFIX}/bin/python
setenv PATH .:${PROTEUS_PREFIX}/bin:${HOME}/bin:${PATH}
setenv LD_LIBRARY_PATH ${PROTEUS_PREFIX}/lib:${TACC_ATLAS_LIB}:/opt/apps/limic2/0.5.4/lib:${LD_LIBRARY_PATH}
setenv MV2_ON_DEMAND_THRESHOLD 2048
