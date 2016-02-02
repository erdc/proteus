#!/bin/bash -x

mpirun -np 2 parun -O petsc.options.asm floating_bar_so.py > logfile.txt
