"""
Module for controlling MPI
"""
import sys
import flcbdfWrappers
initial_communicator = None
isInitialized = 0
def init(petscDatabaseFilename=None,argv=sys.argv):
    global initial_communicator,isInitialized
    if initial_communicator == None:
        #pass
        try:
            import petsc4py
            petsc4py.init(argv)
            isInitialized=1
        except:
            print "WARNING petsc4py import failed!!!"
    if isinstance(petscDatabaseFilename,str):
        comm = flcbdfWrappers.DaetkPetscSys(isInitialized,argv,petscDatabaseFilename)
        isInitialized=1
    else:
        comm = flcbdfWrappers.DaetkPetscSys(isInitialized,argv)
        isInitialized=1
    if initial_communicator == None:
        initial_communicator = comm
    return comm

def get():
    assert initial_communicator != None, "Call Comm.init() once before calling Comm.get()"
    return initial_communicator
