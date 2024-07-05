#! /usr/bin/env python
from glob import glob

pFiles = glob('*_p.py')
for p in pFiles:
    pf = open(p,'r')
    lines = pf.readlines()
    pf.close()
    pf = open(p,'w')
    for line in lines:
        if line.find('PyadhTools') != -1:
            pass
        elif line.find('default_p') != -1:
            line = 'from proteus import *\nfrom proteus.default_p import *\n'
            pf.write(line)
        elif line.find('from VectorTransport') != -1:
            pass
        elif line.find('import AnalyticalSolutions') != -1:
            pass
        else:
            pf.write(line)
pf.close()