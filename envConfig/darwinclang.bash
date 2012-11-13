#basic darwin with default apple compilers (clang/clang++) and mac hpc gfortran #in /usr/local/bin
export PROTEUS_ARCH=darwinclang
export PROTEUS_PREFIX=${PROTEUS}/${PROTEUS_ARCH}
export PROTEUS_PYTHON=${PROTEUS_PREFIX}/Library/Frameworks/Python.framework/Versions/Current/bin/python
export CC=clang
export CXX=clang++
export F90=gfortran
export F77=gfortran
export FC =gfortran
export LD_LIBRARY_PATH=${PROTEUS_PREFIX}/lib
export DYLD_LIBRARY_PATH=${PROTEUS_PREFIX}/lib
export PATH=${PROTEUS_PREFIX}/Library/Frameworks/Python.framework/Versions/Current/bin:${PROTEUS_PREFIX}/bin:${PATH}
