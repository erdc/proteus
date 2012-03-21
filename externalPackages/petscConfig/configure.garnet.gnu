${PROTEUS_PYTHON} ./config/configure.py \
--useThreads=0 \
--prefix=${PROTEUS_PREFIX} \
--PETSC_ARCH=${PROTEUS_ARCH} \
--PETSC_DIR=${PROTEUS}/externalPackages/petsc \
--with-clanguage=C \
--with-cc=cc \
--with-cxx=CC \
--with-fortran=ftn \
--with-pic=1 \
--with-blas-lapack-lib=[/opt/xt-libsci/10.4.9/gnu/lib/45] \
--download-cmake=1 \
--download-metis=1 \
--download-superlu_dist=1 \
--download-parmetis=1 

#--with-shared-libraries=0 \
#--with-fc="ftn" \
#--useThreads=0 \
#\
#--download-blacs=1 \
#--download-scalapack=1 \
#--download-mumps=1 \
#--download-superlu=1 \
#--download-superlu_dist=1 \
#--download-hypre=1
