#basic darwin
setenv MACOSX_DEPLOYMENT_TARGET 10.6
setenv PROTEUS_ARCH darwinhpc
setenv PROTEUS_PREFIX ${PROTEUS}/${PROTEUS_ARCH}
setenv PROTEUS_PYTHON ${PROTEUS_PREFIX}/Library/Frameworks/Python.framework/Versions/Current/bin/python
setenv LD_LIBRARY_PATH ${PROTEUS_PREFIX}/lib
setenv DYLD_LIBRARY_PATH ${PROTEUS_PREFIX}/lib
setenv PATH ${PROTEUS_PREFIX}/Library/Frameworks/Python.framework/Versions/Current/bin:${PROTEUS_PREFIX}/bin:/usr/local/bin:${PATH}
setenv CC  /usr/local/bin/gcc
setenv CXX /usr/local/bin/g++
setenv FC  /usr/local/bin/gfortran
setenv F77 /usr/local/bin/gfortran
setenv F90 /usr/local/bin/gfortran
