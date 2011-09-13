#mac ports
setenv PROTEUS ${HOME}/src/proteus
setenv PROTEUS_ARCH darwinports
setenv PROTEUS_PREFIX ${PROTEUS}/${PROTEUS_ARCH}
setenv PROTEUS_PYTHON ${PROTEUS_PREFIX}/Python.framework/Versions/Current/bin/python
setenv DAETK_DIR ${PROTEUS}/externalPackages/daetk
setenv DAETK_ARCH darwinports
setenv PETSC_DIR ${PROTEUS}/externalPackages/petsc-3.1-p1
setenv PETSC_ARCH darwinports
setenv CC /opt/local/bin/openmpicc
setenv CXX /opt/local/bin/openmpicxx
setenv FC /opt/local/bin/openmpi90
setenv F90 /opt/local/bin/openmpif90
setenv F77 /opt/local/bin/openmpif77
setenv LD_LIBRARY_PATH ${PROTEUS_PREFIX}/lib
setenv DYLD_LIBRARY_PATH ${PROTEUS_PREFIX}/lib
setenv MACOSX_DEPLOYMENT_TARGET 10.6
setenv PATH ${PROTEUS_PREFIX}/Python.framework/Versions/Current/bin:${PROTEUS_PREFIX}/bin:/opt/local/bin:${PATH}
