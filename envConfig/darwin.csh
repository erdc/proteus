#basic darwin
setenv MACOSX_DEPLOYMENT_TARGET 10.6
setenv PROTEUS_ARCH darwin
setenv PROTEUS_PREFIX ${PROTEUS}/${PROTEUS_ARCH}
setenv PROTEUS_PYTHON ${PROTEUS_PREFIX}/Python.framework/Versions/Current/bin/python
setenv F90 gfortran
setenv F77 gfortran
setenv FC  gfortran
setenv LD_LIBRARY_PATH ${PROTEUS_PREFIX}/lib
setenv DYLD_LIBRARY_PATH ${PROTEUS_PREFIX}/lib
setenv PATH ${PROTEUS_PREFIX}/Library/Frameworks/Python.framework/Versions/Current/bin:${PROTEUS_PREFIX}/bin:${PATH}
