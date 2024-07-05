#!/usr/bin/env python
from glob import glob
from os import system
from sys import stdout,argv
oldPrefix = argv[1]
offSet = len(oldPrefix)
newPrefix = argv[2]
pFiles = glob(oldPrefix+'*_p.py')
for p in pFiles:
    nameEnd=p.find('_p.py')
    oldName = p[:nameEnd]
    newName = newPrefix+p[offSet:nameEnd]
    new_p = newName+'_p.py'
    system('svn mv '+p+' '+new_p)
    #print 'svn mv '+p+' '+new_p
    nFiles = glob(p[:nameEnd]+'*_n.py')
    for n in nFiles:
        nameEnd_n = n.find('_n.py')
        newName_n = newPrefix+n[offSet:nameEnd_n]
        new_n = newName_n+'_n.py'
        #print 'svn mv '+n+' '+new_n
        system('svn mv '+n+' '+new_n)
        nf = open(new_n,'r')
        #nf = open(n,'r')
        lines = nf.readlines()
        nf.close()
        nf = open(new_n,'w')
        for line in lines:
            if line.find('from') != -1:
                if line.find('_p') != -1:
                    line = 'from '+newName+'_p import *\n'
            #stdout.write(line)
            nf.write(line)
nf.close()