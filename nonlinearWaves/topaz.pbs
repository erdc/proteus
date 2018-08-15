#!/bin/bash
#PBS -A ERDCV00898R40
#PBS -l walltime=001:00:00
#PBS -l select=1:ncpus=36:mpiprocs=36
#PBS -q debug
#PBS -N nlwaves
#PBS -j oe
#PBS -l application=proteus
#PBS -V
#PBS -m eba
#PBS -M cekees@gmail.com
cd $PBS_O_WORKDIR
mkdir $WORKDIR/$PBS_JOBNAME.$PBS_JOBID
cp *.py $WORKDIR/$PBS_JOBNAME.$PBS_JOBID
cp topaz.pbs $WORKDIR/$PBS_JOBNAME.$PBS_JOBID
cp $HOME/air-water-vv/inputTemplates/petsc.options.superlu_dist $WORKDIR/$PBS_JOBNAME.$PBS_JOBID
cd  $WORKDIR/$PBS_JOBNAME.$PBS_JOBID
mpiexec_mpt -n ${BC_MPI_TASKS_ALLOC}  parun nonlinear_waves_so.py -p -l 5 -v -O petsc.options.superlu_dist
