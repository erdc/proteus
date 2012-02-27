${PROTEUS_PYTHON} ./config/configure.py \
--with-pic \
--with-clanguage=C \
--with-cc=mpicc \
--with-cxx=mpicxx \
--with-fc=mpif90 \
--with-mpi-compilers \
--with-shared-libraries \
--download-spooles=1 \
--download-superlu_dist=1 \
--download-hypre=1 \
--with-blas-lapack-dir=${TACC_ATLAS_LIB} \
--download-metis=1 \
--download-parmetis=1 \
--download-cmake=1 \
--prefix=${PROTEUS_PREFIX} \
--PETSC_ARCH=${PROTEUS_ARCH} \
--PETSC_DIR=${PROTEUS}/externalPackages/petsc
#--with-metis=1 \
#--with-metis-dir=/corral/hurricane/cekees/src-lonestar/proteus/externalPackages/petsc/lonestar.gnu \
#--with-parmetis-dir=/corral/hurricane/cekees/src-lonestar/proteus/externalPackages/petsc/lonestar.gnu \



