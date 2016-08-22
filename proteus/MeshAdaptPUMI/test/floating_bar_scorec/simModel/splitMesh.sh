#!/bin/bash -x
#$1 is number of processors
#$2 is the model name
mpirun -np $1 /lore/zhanga/core-sim/build/test/zsplit $2.dmg $2.smb $1-Procs/ $1
