#!/bin/bash +x

mpirun -np 2 parun -O petsc.options splashcube_so.py > logfile.txt
