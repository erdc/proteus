#!/usr/bin/env python
import shelve
from heavyTest import *

lightData = shelve.open('lightData',flag='c',protocol=0)
heavyData = shelve.open('heavyData',flag='c',protocol=0)

for sol in lightData.values():
    print "T=",sol.T
    sol.reinitialize(heavyData)
    print "sol",sol.sol
