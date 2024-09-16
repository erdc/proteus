#import sys
#from proteus import Comm
import os
import sys
import socket
import pickle
import numpy
import proteus
#remove blankent import statements until after Comm initialized for petsc4py
# from proteus import *
from proteus import Profiling,Comm
from warnings import *
import optparse
import sys
import pstats
import pdb
import petsc4py
Profiling.openLog("proteus.log",7)
print(sys.argv[:1])
Comm.argv = sys.argv[:1]
comm = Comm.init()
comm = Comm.init()
#petsc4py.init(sys.argv[:1])
print(comm.rank(),comm.size())
print("Hello World from",comm.rank())
from proteus import *
Profiling.procID=comm.rank()
