""" basic utilities in Proteus
"""

import os

def get_include_dir():
    return os.path.dirname(os.path.realpath(__file__))

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