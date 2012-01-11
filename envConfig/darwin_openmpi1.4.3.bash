#matt's version of openmpicc
export PROTEUS_ARCH=darwin_openmpi1.4.3
export PROTEUS_PREFIX=${PROTEUS}/${PROTEUS_ARCH}
export PROTEUS_PYTHON=${PROTEUS_PREFIX}/Library/Frameworks/Python.framework/Versions/Current/bin/python
export CC=/opt/openMPI-1.4.3/bin/mpicc
export CXX=/opt/openMPI-1.4.3/bin/mpicxx
export FC=/opt/openMPI-1.4.3/bin/mpi90
export F90=/opt/openMPI-1.4.3/bin/mpif90
export F77=/opt/openMPI-1.4.3/bin/mpif77
export LD_LIBRARY_PATH=${PROTEUS_PREFIX}/lib
export DYLD_LIBRARY_PATH=${PROTEUS_PREFIX}/lib
export MACOSX_DEPLOYMENT_TARGET=10.6
export PATH=${PROTEUS_PREFIX}/Library/Frameworks/Python.framework/Versions/Current/bin:${PROTEUS_PREFIX}/bin:/opt/openMPI-1.4.3/bin:${PATH}
