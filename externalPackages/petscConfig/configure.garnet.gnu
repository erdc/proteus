${PROTEUS_PYTHON} ./config/configure.py \
--useThreads=0 \
--with-make-np=1 \
--prefix=${PROTEUS_PREFIX} \
--PETSC_ARCH=${PROTEUS_ARCH} \
--PETSC_DIR=${PROTEUS}/externalPackages/petsc \
--with-clanguage=C \
--with-mpi-compilers=1 \
--CC=cc \
--CXX=CC \
--with-fortran=0 \
--with-pic=1 \
--with-blas-lapack-dir=[/opt/xt-libsci/10.4.9/gnu/lib/45] \
--with-blas-lapack-lib=[-lsci] \
--download-cmake=1 \
--download-metis=1 \
--download-parmetis=1 

#--with-pthread=0 \
#--download-superlu_dist=1 \
#--FC=ftn \


#--sharedLibraryFlags="-dynamic" \

#--with-fortran=ftn \
#--useThreads=0 \
#\
#--download-blacs=1 \
#--download-scalapack=1 \
#--download-mumps=1 \
#--download-superlu=1 \
#--download-superlu_dist=1 \
#--download-hypre=1
