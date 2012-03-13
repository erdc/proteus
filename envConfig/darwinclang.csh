#basic darwin
setenv MACOSX_DEPLOYMENT_TARGET 10.7
setenv PROTEUS_ARCH darwinclang
setenv PROTEUS_PREFIX ${PROTEUS}/${PROTEUS_ARCH}
setenv PROTEUS_PYTHON ${PROTEUS_PREFIX}/Library/Frameworks/Python.framework/Versions/Current/bin/python
setenv CC clang
setenv CXX clang++
setenv F90 gfortran
setenv F77 gfortran
setenv FC  gfortran
setenv LD_LIBRARY_PATH ${PROTEUS_PREFIX}/lib:${PROTEUS_PREFIX}/lib/paraview-3.14
setenv DYLD_LIBRARY_PATH ${PROTEUS_PREFIX}/lib:${PROTEUS_PREFIX}/lib/paraview-3.14
setenv PATH ${PROTEUS_PREFIX}/Library/Frameworks/Python.framework/Versions/Current/bin:${PROTEUS_PREFIX}/bin:${PATH}
