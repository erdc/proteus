from proteus.default_so import *
import dambreak

from proteus.SplitOperator import Sequential_FixedStep_Simple, defaultSystem

pnList = [("clsvof_p",               "clsvof_n"),#0
          ("twp_navier_stokes_p", "twp_navier_stokes_n"),#1
          ("pressureincrement_p", "pressureincrement_n"),#2
          ("pressure_p", "pressure_n"),#3
          ("pressureInitial_p", "pressureInitial_n")]#4
dambreak.VOS_model=None
dambreak.SED_model=None
dambreak.VOF_model=None
dambreak.LS_model=None
dambreak.RD_model=None
dambreak.MCORR_model=None
dambreak.CLSVOF_model=0
dambreak.V_model=1
dambreak.PINC_model=2
dambreak.PRESSURE_model=3
dambreak.PINIT_model=4

name = "dambreak"

modelSpinUpList = [dambreak.PINIT_model]

class Sequential_MinAdaptiveModelStepPS(Sequential_MinAdaptiveModelStep):
    def __init__(self,modelList,system=defaultSystem,stepExact=True):
        Sequential_MinAdaptiveModelStep.__init__(self,modelList,system,stepExact)
        self.modelList = modelList[:len(pnList)-1]


systemStepControllerType = Sequential_MinAdaptiveModelStepPS

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,dambreak.dt_init]+[i*dambreak.dt_fixed for i in range(1,dambreak.nDTout+1)]

info = open("TimeList.txt","w")
#archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
