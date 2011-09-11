#basic darwin64
setenv PROTEUS ${HOME}/src/proteus
setenv PROTEUS_ARCH darwin64
setenv PROTEUS_PREFIX ${PROTEUS}/${PROTEUS_ARCH}
setenv PROTEUS_PYTHON ${PROTEUS_PREFIX}/Python.framework/Versions/Current/bin/python
setenv DAETK_DIR ${PROTEUS}/externalPackages/daetk
setenv DAETK_ARCH darwin64
setenv PETSC_DIR ${PROTEUS}/externalPackages/petsc-3.1-p1
setenv PETSC_ARCH darwin64
setenv CC /usr/bin/mpicc
setenv CFLAGS "-arch x86_64"
setenv CPPFLAGS "-arch x86_64"
setenv CXX /usr/bin/mpicxx
setenv CXXFLAGS "-arch x86_64"
setenv FC  /opt/local/bin/openmpif90
setenv FFLAGS "-m64"
setenv F77 /opt/local/bin/openmpif77
setenv F90 /opt/local/bin/openmpif90
setenv LD_LIBRARY_PATH ${PROTEUS_PREFIX}/lib
setenv DYLD_LIBRARY_PATH ${PROTEUS_PREFIX}/lib
setenv MACOSX_DEPLOYMENT_TARGET 10.6
setenv PATH ${PROTEUS_PREFIX}/Python.framework/Versions/Current/bin:${PROTEUS_PREFIX}/bin:${PATH}
