from proteus.default_so import *
import NS_convergence

from proteus.SplitOperator import Sequential_FixedStep_Simple, defaultSystem

pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),#0
          ("pressure_p", "pressure_n")] #1
NS_convergence.VOF_model=None
NS_convergence.VOS_model=None
NS_convergence.SED_model=None
NS_convergence.V_model=0
NS_convergence.PINC_model=None
NS_convergence.PRESSURE_model=1
NS_convergence.PINIT_model=None

name = "NS_convergence"

class Sequential_MinAdaptiveModelStepPS(Sequential_MinAdaptiveModelStep):
    def __init__(self,modelList,system=defaultSystem,stepExact=True):
        Sequential_MinAdaptiveModelStep.__init__(self,modelList,system,stepExact)
        self.modelList = modelList[:len(pnList)-1]
        
systemStepControllerType = Sequential_MinAdaptiveModelStepPS
#systemStepControllerType = Sequential_FixedStep
#dt_system_fixed = 0.5*NS_convergence.he

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,NS_convergence.dt_init]+[i*NS_convergence.dt_fixed for i in range(1,NS_convergence.nDTout+1)]
