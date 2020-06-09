#!/bin/bash 

#generating the PUMI mesh needs to be run in serial (1-Processor / core)
#parun -l 5 --TwoPhaseFlow --genPUMI dambreak.py -v -D "R00_PUMI"
aprun -n 1 parun -l 5 --TwoPhaseFlow --genPUMI dambreak.py -v -D "R00_PUMI"

