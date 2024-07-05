#! /usr/bin/env python
import sys
import os
import glob

pFiles = glob.glob('*_p.py')
caseDict = {}
for pf in pFiles:
    caseDict[pf] = set(glob.glob(pf[:-5]+'*_n.py'))
#fix cases were problem name is a subset of some other problem name
for pf1 in pFiles:
    for pf2 in pFiles:
        if pf2.find(pf1[:-4]):
            nf1Set=set(glob.glob(pf1[:-5]+'*_n.py'))
            caseDict[pf2] -= nf1Set
for pf in pFiles:
    print(pf)
    print(caseDict[pf])
for p,nList in caseDict.items():
    if len(nList) == 0:
        sys.stdout.write("\n----------------Skipping "+p+".  No n file----------------------\n")
        sys.stdout.flush()
    else:
        for n in nList:
            args = ('proteusRun.py',p,n,'-l 4','-b','runAllBatch.py')
            sys.stdout.write("\n----------------Running  "+p+"---"+n+"\n")
            sys.stdout.flush()
            os.spawnvpe(os.P_WAIT,'proteusRun.py',args,os.environ)
