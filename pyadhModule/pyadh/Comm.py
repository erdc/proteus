"""
Module for controlling MPI
TODO:
  access petsc4py PetscSys instead of through flcbdfWrappers
  pass petscDatabaseFilename to petsc4py
"""
import sys
import flcbdfWrappers
#mwf debug
initial_communicator = None
def init(petscDatabaseFilename=None,argv=sys.argv):
    try:
        import petsc4py
        petsc4py.init(argv)
                
    except:
        print "WARNING petsc4py import failed!!!" 
    global initial_communicator
    if isinstance(petscDatabaseFilename,str):
        comm = flcbdfWrappers.DaetkPetscSys(argv,petscDatabaseFilename)
    else:
        comm = flcbdfWrappers.DaetkPetscSys(argv)
    if initial_communicator == None:
        initial_communicator = comm
    return comm

def get():
    if initial_communicator != None:
        return initial_communicator
    return flcbdfWrappers.getDaetkPetscSys()
