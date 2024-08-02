import proteus.default_so
from importlib import reload
reload(proteus.default_so)
from proteus.default_so import *
from . import NS_hotstart
reload(NS_hotstart)

from proteus.SplitOperator import Sequential_FixedStep_Simple, defaultSystem

pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),#0
          ("pressureincrement_p", "pressureincrement_n"),#1
          ("pressure_p", "pressure_n"),
          ("pressureInitial_p", "pressureInitial_n")]#3

NS_hotstart.VOF_model=None
NS_hotstart.VOS_model=None
NS_hotstart.SED_model=None
NS_hotstart.V_model=0
NS_hotstart.PINC_model=1
NS_hotstart.PRESSURE_model=2
NS_hotstart.PINIT_model=3

name = "NS_hotstart"

# class Sequential_MinAdaptiveModelStepPS(Sequential_MinAdaptiveModelStep):
#     def __init__(self,modelList,system=defaultSystem,stepExact=True):
#         Sequential_MinAdaptiveModelStep.__init__(self,modelList,system,stepExact)
#         self.modelList = modelList[:len(pnList)-1]
#         
# systemStepControllerType = Sequential_MinAdaptiveModelStepPS

modelSpinUpList = [NS_hotstart.PINIT_model]

class Sequential_MinAdaptiveModelStepPS(Sequential_FixedStep):
    def __init__(self,modelList,system=defaultSystem,stepExact=True):
        Sequential_FixedStep.__init__(self,modelList,system,stepExact)
        self.modelList = modelList[:len(pnList)-1]

dt_system_fixed = NS_hotstart.dt_fixed
systemStepControllerType = Sequential_MinAdaptiveModelStepPS

needEBQ_GLOBAL = False
needEBQ = False

tnList = NS_hotstart.tnList
