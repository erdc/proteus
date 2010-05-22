from pyadh import Comm
import sys
comm = Comm.init()

#this wasn't working until I let pyadh handle the actual MPI init and finalize
from pyadh import ADH
#het_column0 = ADH.ADH_OneLevelTransport(adhInput="het_column_8",runname="test1")
het_column0 = ADH.ADH_OneLevelTransport(adhInput="angle_bt",runname="test1")
#het_column0.calculateSolution("stringNotYetUsed")


