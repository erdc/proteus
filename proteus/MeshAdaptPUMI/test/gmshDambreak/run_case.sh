#!/bin/bash +x

#valgrind --log-file=valgrind%p --leak-check=full parun couette_so.py -l5 -v
#parun dambreak_so.py -v -D runTest_1
mpirun -np 4 parun -O petsc.options dambreak_so.py -v -D runTest2_parallel
