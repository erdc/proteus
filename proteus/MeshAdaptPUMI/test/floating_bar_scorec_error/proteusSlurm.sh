#!/bin/bash -x 


#SBATCH -J FloatingBar
#SBATCH -o FloatingBar.%j.log
#SBATCH -e FloatingBar.%j.error
#SBATCH -n 8
#SBATCH -t 2160:00
#SBATCH -p ib
#SBATCH -N 1

NPROCS=8

mpirun -n $NPROCS parun -O petsc.options.asm floating_bar_so.py
