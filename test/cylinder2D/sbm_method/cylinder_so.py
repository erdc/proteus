import proteus.default_so
from importlib import reload
reload(proteus.default_so)
from proteus.default_so import *
try:
    from . import cylinder
except:
    import cylinder
reload(cylinder)

from proteus.SplitOperator import Sequential_FixedStep_Simple, defaultSystem

pnList = cylinder.pnList

if cylinder.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "cylinder"
soname=name
#modelSpinUpList = [cylinder.VOF_model, cylinder.LS_model, cylinder.V_model, cylinder.PINIT_model]
modelSpinUpList = [cylinder.PINIT_model]

# class Sequential_MinAdaptiveModelStepPS(Sequential_MinAdaptiveModelStep):
#     def __init__(self,modelList,system=defaultSystem,stepExact=True):
#         Sequential_MinAdaptiveModelStep.__init__(self,modelList,system,stepExact)
#         self.modelList = modelList[:len(pnList)-1]


class Sequential_MinAdaptiveModelStepPS(Sequential_FixedStep):
    def __init__(self,modelList,system=defaultSystem,stepExact=True):
        Sequential_FixedStep.__init__(self,modelList,system,stepExact)
        self.modelList = modelList[:len(pnList)-1]

dt_system_fixed = cylinder.dt_fixed
systemStepControllerType = Sequential_MinAdaptiveModelStepPS

# systemStepControllerType = Sequential_FixedStep #Sequential_FixedStep_Simple # uses time steps in so.tnList
# dt_system_fixed = 0.01; 
systemStepExact=False


needEBQ_GLOBAL = False
needEBQ = False

# tnList = [0.0,cylinder.dt_init]+[i*cylinder.dt_fixed for i in range(1,cylinder.nDTout+1)]
tnList = cylinder.tnList
info = open("TimeList.txt","w")
#archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
