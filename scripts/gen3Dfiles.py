#!/usr/bin/env python
from glob import glob
from os import system
from sys import stdout,argv
problem2d = argv[1]
problem3d = problem2d[:-2]+'3d'
pnFiles = glob('*'+problem2d+'*.py')
print('*'+problem2d+'.py',pnFiles)
for f in pnFiles:
    nameEnd=f.find(problem2d)
    newName= f[:nameEnd]+problem3d+f[nameEnd+len(problem2d):]
    cmd = 'svn cp '+f+' '+newName
    system(cmd)
