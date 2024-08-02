#!/usr/bin/env python
from glob import glob
from sys import stdout
nFiles = glob('*_n.py')
for n in nFiles:
    nf = open(n,'r')
    lines = nf.readlines()
    nf.close()
    nf = open(n,'w')
    for line in lines:
        start=line.find('2d_p import *')
        if start != -1:
            line = line[:start+2]+'_p import *\n'
        start = line.find('1d_p import *')
        if start != -1:
            line = line[:start+2]+'_p import *\n'
        start = line.find('ladr')
        if start != -1:
            line = line[:start]+'ladr'+line[start+6:]+'\n'
        #stdout.write(line)
        nf.write(line)
    nf.close()