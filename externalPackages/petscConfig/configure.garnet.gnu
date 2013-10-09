#!/bin/bash
unset CC
unset CXX
unset CXX
${PROTEUS_PYTHON} ./config/configure.py \
CC=cc \
CXX=CC \
FC=ftn \
--with-debugging=0 \
--useThreads=0 \
--prefix=${PROTEUS_PREFIX} \
--PETSC_ARCH=${PROTEUS_ARCH} \
--PETSC_DIR=${PROTEUS}/externalPackages/petsc-dev \
--with-clanguage=C \
--with-mpi-compilers=1 \
--with-blas-lapack-lib="-lsci_gnu" \
--with-pic=1 \
--download-cmake=1 \
--download-metis=1 \
--download-superlu=1 \
--download-superlu_dist=1 \
--download-parmetis=1 \
--download-blacs=1 \
--download-scalapack=1 \
--download-mumps=1 \
--download-hypre=1
