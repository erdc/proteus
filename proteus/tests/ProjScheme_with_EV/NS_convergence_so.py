from proteus.default_so import *
from . import NS_convergence

from proteus.SplitOperator import Sequential_FixedStep_Simple, defaultSystem

pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),#0
          ("pressureincrement_p", "pressureincrement_n"),#1
          ("pressure_p", "pressure_n"),#2
          ("pressureInitial_p", "pressureInitial_n")]#3

NS_convergence.VOF_model=None
NS_convergence.VOS_model=None
NS_convergence.SED_model=None
NS_convergence.V_model=0
NS_convergence.PINC_model=1
NS_convergence.PRESSURE_model=2
NS_convergence.PINIT_model=3

name = "NS_convergence"

class Sequential_MinAdaptiveModelStepPS(Sequential_MinAdaptiveModelStep):
    def __init__(self,modelList,system=defaultSystem,stepExact=True):
        Sequential_MinAdaptiveModelStep.__init__(self,modelList,system,stepExact)
        self.modelList = modelList[:len(pnList)-1]
        
systemStepControllerType = Sequential_MinAdaptiveModelStepPS

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,NS_convergence.dt_init]+[i*NS_convergence.dt_fixed for i in range(1,NS_convergence.nDTout+1)]
