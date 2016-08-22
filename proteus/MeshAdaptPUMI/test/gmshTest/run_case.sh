#!/bin/bash +x

#valgrind --log-file=valgrind%p --leak-check=full parun couette_so.py -l5 -v
parun couette_so.py -l5 -v > logfile.txt
