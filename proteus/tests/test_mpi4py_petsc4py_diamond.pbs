#! /bin/bash
#PBS -A ERDCV00898ENQ
#PBS -l walltime=001:00:00
#PBS -N test
#PBS -q debug
#PBS -j oe
#PBS -l select=2:ncpus=8:mpiprocs=8
#PBS -l place=scatter:excl
#PBS -l application=proteus
cd ${PROTEUS}/externalPackages/mpi4py/test
mpirun -np 16 ${PROTEUS_PYTHON} runtests.py --exclude test_spawn
cd ${PROTEUS}/externalPackages/petsc4py/test
mpirun -np 16 ${PROTEUS_PYTHON} runtests.py --exclude test_da --exclude test_obj --exclude test_vec --exclude test_ts -v
echo "done testing mpi4py and petsc4py"
