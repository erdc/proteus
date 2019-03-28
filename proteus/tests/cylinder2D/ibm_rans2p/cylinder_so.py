from proteus import *
from proteus.default_so import *
import cylinder

from proteus.SplitOperator import Sequential_FixedStep_Simple, defaultSystem


pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),#0
          ]#3


name = "cylinder"

# class Sequential_MinAdaptiveModelStepPS(Sequential_MinAdaptiveModelStep):
#     def __init__(self,modelList,system=defaultSystem,stepExact=True):
#         Sequential_MinAdaptiveModelStep.__init__(self,modelList,system,stepExact)
#         self.modelList = modelList[:len(pnList)-1]


class Sequential_MinAdaptiveModelStepPS(Sequential_FixedStep):
    def __init__(self,modelList,system=defaultSystem,stepExact=True):
        Sequential_FixedStep.__init__(self,modelList,system,stepExact)
        self.modelList = modelList[:len(pnList)]

dt_system_fixed = cylinder.dt_fixed
systemStepControllerType = Sequential_MinAdaptiveModelStepPS

# stepController = StepControl.Min_dt_cfl_controller
systemStepExact=False

# systemStepControllerType = Sequential_FixedStep #Sequential_FixedStep_Simple # uses time steps in so.tnList
# dt_system_fixed = 0.01; 
# systemStepExact=False;

tnList = cylinder.tnList

info = open("TimeList.txt","w")
#archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
