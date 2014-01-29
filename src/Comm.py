"""
Module for controlling MPI
"""
import sys
from proteus import flcbdfWrappers
import petsc4py

initial_communicator = None
isInitialized = 0
def set_isInitialized():
    global isInitialized
    isInitialized=1

def init(petscDatabaseFilename=None,argv=sys.argv):
    global initial_communicator,isInitialized
    if isInitialized is 0:
        petsc4py.init(argv)
        isInitialized=1
    else:
        from petsc4py import PETSc
        narg = len(argv)
        for i,s in enumerate(argv):
            if len(s) > 0 and s[0] is '-':
                name = s
                if i+1 < narg and  argv[i+1][0] is not '-':
                    value = argv[i+1]
                    PETSc.Options.SetValue(name,value)
    if isinstance(petscDatabaseFilename,str):
        comm = flcbdfWrappers.DaetkPetscSys(isInitialized,argv,petscDatabaseFilename)
    else:
        comm = flcbdfWrappers.DaetkPetscSys(isInitialized,argv)
    if initial_communicator == None:
        initial_communicator = comm
    return comm

def get():
    assert initial_communicator != None, "Call Comm.init() once before calling Comm.get()"
    return initial_communicator
