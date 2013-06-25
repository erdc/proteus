import numpy
from proteus import Comm
comm=Comm.init()
if comm.isMaster():
    numpy.test(verbose=2)
    print("PEXPECT_EXIT")
