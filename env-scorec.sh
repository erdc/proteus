module load mpich3/3.1.2-thread-multiple git python cmake simmetrix/simModSuite gcc/4.7.0
export HDPATH=/lore/zhanga/HDPATH2
export PROTEUS=$PWD
export PROTEUS_ARCH=linux2-scorec
export PROTEUS_PREFIX=$PROTEUS/$PROTEUS_ARCH
export PYTHONPATH=$PROTEUS_PREFIX/lib/python2.7/site-packages
export PROTEUS_PYTHON=$PROTEUS_PREFIX/bin/python
export PATH=$PROTEUS_PREFIX/bin:$PATH
export MPI_ROOT=/usr/local/mpich3/3.1.2-thread-multiple
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PROTEUS/linux2-scorec/lib/:$MPI_ROOT/lib
export EDITOR=vim
export PATH=$PROTEUS/hashdist/bin:$PATH

#GMSH support

export PATH=/users/zhanga/gmsh-2.10.1-Linux/bin:$PATH

#export PETSC_DIR=/users/zhanga/Projects/PETSc/petsc
#export PETSC_ARCH=arch-linux2-c-debug

#export PETSC_INCLUDE_DIR=$PETSC_DIR/$PETSC_ARCH/include
#export PETSC_LIB_DIR=$PETSC_DIR/$PETSC_ARCH/lib
