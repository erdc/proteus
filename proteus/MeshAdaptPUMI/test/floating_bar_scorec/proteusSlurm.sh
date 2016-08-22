#!/bin/bash -x 

HOSTFILE=/tmp/hosts.$SLURM_JOB_ID
srun hostname -s > $HOSTFILE

export SIM_LICENSE_FILE=/fasttmp/chitak/license.txt

#SBATCH -J FloatingBar
#SBATCH -o FloatingBar.%j.log
#SBATCH -e FloatingBar.%j.error
#SBATCH -n 8
#SBATCH -t 2160:00
#SBATCH -p ib
#SBATCH -N 1

if [ -z "$SLURM_NPROCS" ] ; then
  if [ -z "$SLURM_NTASKS_PER_NODE" ] ; then
    SLURM_NTASKS_PER_NODE=1
  fi
  SLURM_NPROCS=$(( $SLURM_JOB_NUM_NODES * $SLURM_NTASKS_PER_NODE ))
fi

/usr/local/mpich3/3.1.2-thread-multiple/bin/mpirun -machinefile $HOSTFILE -np $SLURM_NPROCS parun -O petsc.options.asm floating_bar_so.py

rm /tmp/hosts.$SLURM_JOB_ID
