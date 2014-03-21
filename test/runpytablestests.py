import tables
from proteus import Comm
comm = Comm.init()
if comm.isMaster():
    tables.test(verbose=2)
    print("PEXPECT_EXIT")
