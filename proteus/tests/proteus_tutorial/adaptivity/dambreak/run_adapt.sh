#!/bin/bash 

#mpirun -np 4 parun -l 5 --TwoPhaseFlow dambreak_PUMI.py -v -D "R00_PUMI"
aprun -n 4 parun -l 5 --TwoPhaseFlow dambreak_PUMI.py -v -D "R00_PUMI"

