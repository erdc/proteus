"""
Module for controlling MPI

.. inheritance-diagram:: proteus.Comm
   :parts: 1
"""
from __future__ import absolute_import

from builtins import range
from builtins import object
import ctypes
import sys

import petsc4py
from .Profiling import logEvent


# Special workaround for broken MPI on certain Cray systems
from . import config
mpi_preload_libs=[]
for lib in config.PROTEUS_PRELOAD_LIBS:
    mpi_preload_libs.append(ctypes.CDLL(lib,mode=ctypes.RTLD_GLOBAL))

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
        if not OptDB.hasName("options_left"):
            OptDB.setValue("options_left",False)
    petsc4py.PETSc.Log.Stage('proteus').push()
    new_comm = Comm()
    if not comm:
        comm = new_comm
    return new_comm

def get():
    if comm is None:
        logEvent("Comm.get called before init, init is being called for you.", 3)
        return init()
    return comm

class Comm(object):
    """Proteus wrapper around PETSc/MPI communicators

    This is a very simple class that provides compatibility with older Proteus comm objects.
    """

    def __init__(self):
        from petsc4py import PETSc
        self.comm = PETSc.COMM_WORLD
        self.mpi4py_comm = PETSc.COMM_WORLD.tompi4py()
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

    def globalSum(self, value):
        from mpi4py import MPI
        rvalue = self.mpi4py_comm.allreduce(sendobj = value,
                                            op = MPI.SUM)
        return rvalue
    
    def globalMax(self, value):
        from mpi4py import MPI
        rvalue = self.mpi4py_comm.allreduce(sendobj = value,
                                            op = MPI.MAX)
        return rvalue

    def globalMin(self, value):
        from mpi4py import MPI
        rvalue = self.mpi4py_comm.allreduce(sendobj = value,
                                            op = MPI.MIN)
        return rvalue
    
def globalSum(value):
    return comm.globalSum(value)

def globalMax(value):
    return comm.globalMax(value)

def globalMin(value):
    return comm.globalMin(value)
