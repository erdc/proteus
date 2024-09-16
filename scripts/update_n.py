#! /usr/bin/env python
from glob import glob

nFiles = glob('*_n.py')

for n in nFiles:
    nf = open(n,'r')
    lines = nf.readlines()
    nf.close()
    nf = open(n,'w')
    for line in lines:
        if line.find('from PyadhTools import *') != -1:
            pass
        elif line.find('from default_n import *') != -1:
            line = 'from proteus import *\nfrom proteus.default_n import *\n'
            nf.write(line)
        elif line.find('conservativeFlux = \'pwl\'') != -1:
            line = 'conservativeFlux = {0:\'pwl\'}'
            nf.write(line)
        elif line.find('shockCapturing = {') != -1:
            line = 'shockCapturing = None\n'
            nf.write(line)
        elif line.find('shockCapturingDiffusion = {') != -1:
            pass
        elif line.find('from defaultHJ') != -1:
            pass
        elif line.find('subgridError = HJstab') != -1:
            line = 'subgridError = HamiltonJacobi_ASGS(coefficients,nd)\n'
            nf.write(line)
        elif line.find('subgridError = HamiltonJacobi_ASGS') != -1:
            line = 'subgridError = HamiltonJacobi_ASGS(coefficients,nd)\n'
            nf.write(line)
        else:
            nf.write(line)
    nf.close()