#!/bin/bash
#PBS -A HPCMO34111002
#PBS -l walltime=024:00:00
#PBS -l select=64:ncpus=32:mpiprocs=32
#PBS -l place=scatter:excl
#PBS -q standard
#PBS -N marin
#PBS -j oe
#PBS -l application=proteus
#PBS -V
#PBS -m eba
##PBS -M myemail@mydomain
#setup modules and check proteus version
#. /opt/cray/pe/modules/default/etc/modules.sh
#export MODULEPATH=/p/app/unsupported/proteus/modulefiles:${MODULEPATH}
#module load proteus/master
#module load proteus/1.7.0
which parun
#create a working directory and copy over the inputs use
cd $PBS_O_WORKDIR
mkdir $WORKDIR/$PBS_JOBNAME.$PBS_JOBID
cp marin.py $WORKDIR/$PBS_JOBNAME.$PBS_JOBID
cp marin.pbs $WORKDIR/$PBS_JOBNAME.$PBS_JOBID
#change into the work directory and run
cd  $WORKDIR/$PBS_JOBNAME.$PBS_JOBID
aprun -n ${BC_MPI_TASKS_ALLOC}  parun -F -l 5 --TwoPhaseFlow marin.py -C "he=0.0125"
