#!/bin/bash -x
HOSTFILE=/tmp/hosts.$SLURM_JOB_ID

srun hostname -s > $HOSTFILE

if [ -z "$SLURM_NPROCS" ] ; then
  if [ -z "$SLURM_NTASKS_PER_NODE" ] ; then
     SLURM_NTASKS_PER_NODE=1
  fi  
	  SLURM_NPROCS=$(( $SLURM_JOB_NUM_NODES * $SLURM_NTASKS_PER_NODE ))
fi
echo 'Starting job'
date
module load python/2.7.6
 cd /fasttmp/chitak/ProtoScorec/proteus
  . env.sh
  cd -
 export SIM_LICENSE_FILE=/fasttmp/chitak/license.txt
 /usr/local/openmpi/latest/bin/mpirun -machinefile $HOSTFILE -np  $SLURM_NPROCS $PWD/parun marin_so.py > log 
 #/usr/local/openmpi/latest/bin/mpirun -machinefile $HOSTFILE -np  $SLURM_NPROCS valgrind --suppressions==mpi.supp
 /fasttmp/chitak/Athena/athenavms.Cmake/athenavms --local --nprocs 64
 --stk -c t4_bar_pumi.cntl --rpi t4_bar_pumi_part.sms --exodus
  out_box.exo --acis newmodel.sat |& tee log_parbox.log

cat /tmp/hosts.$SLURM_JOB_ID
rm /tmp/hosts.$SLURM_JOB_ID
echo 'Job completed'
date
