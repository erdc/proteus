"""
Module for controlling MPI
"""
import sys
import petsc4py

comm = None
argv = sys.argv
petscInitialized = False

def init():
    global comm, petscInitialized
    if not petscInitialized:
        petsc4py.init(argv)
        petscInitialized = True
    else:
        from petsc4py import PETSc
        OptDB=PETSc.Options()
        narg = len(argv)
        for i,s in enumerate(argv):
            if len(s) > 0 and s[0] is '-':
                name = s
                if i+1 < narg and  argv[i+1][0] is not '-':
                    value = argv[i+1]
                    OptDB.setValue(name,value)
    new_comm = Comm()
    if not comm:
        comm = new_comm
    return new_comm

def get():
    assert(comm, "Call Comm.init() once before calling Comm.get()")
    return comm

class Comm():
    """Proteus wrapper around PETSc/MPI communicators

    This is a very simple class that provides compatibility with older Proteus comm objects.
    """

    def __init__(self):
        from petsc4py import PETSc
        self.comm = PETSc.COMM_WORLD

    def isInitialized(self):
        return True

    def rank(self):
        """Return position in this communicator (0-indexed)"""
        return self.comm.rank

    def size(self):
        """Return number of processes on this communicator"""
        return self.comm.size

    def barrier(self):
        """Block all processes on this communicator until all have called this method"""
        return self.comm.Barrier()

    def isMaster(self):
        """Return true if this process is the 'master' for this communicator"""
        return self.comm.rank == 0

    def beginSequential(self, ng=1):
        """Begin a block of code to be completed sequentially by each process
        """

        if ng != 1:
            raise NotImplementedError("Only ng = 1 allowed")

        for i in range(self.comm.rank):
            self.comm.Barrier()
        return

    def endSequential(self, ng=1):
        """End a block of code to be completed sequentially by each process
        """

        if ng != 1:
            raise NotImplementedError("Only ng = 1 allowed")

        for i in range(self.comm.size - self.comm.rank - 1):
            self.comm.Barrier()
        return