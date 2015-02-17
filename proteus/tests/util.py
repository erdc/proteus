import os
import sys
from proteus import Comm, Profiling

def setup_profiling():
    comm = Comm.get()
    Profiling.procID = comm.rank()
    Profiling.logLevel = 10
    Profiling.logFile = sys.stdout
    Profiling.logAllProcesses = True

def silent_rm(filename):
    try:
        os.remove(filename)
    except OSError:
        pass