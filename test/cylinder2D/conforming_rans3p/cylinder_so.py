import proteus.default_so
from importlib import reload
reload(proteus.default_so)
from proteus.default_so import *
try:
    from . import cylinder3p
except:
    import cylinder3p

from proteus.SplitOperator import Sequential_FixedStep_Simple, defaultSystem

pnList = cylinder3p.pnList

# if cylinder3p.useRANS > 0:
#     pnList.append(("kappa_p",
#                    "kappa_n"))
#     pnList.append(("dissipation_p",
#                    "dissipation_n"))
name = "cylinder"
soname=name
#modelSpinUpList = [cylinder3p.VOF_model, cylinder3p.LS_model, cylinder3p.V_model, cylinder3p.PINIT_model]
modelSpinUpList = [cylinder3p.PINIT_model]

class Sequential_MinAdaptiveModelStepPS(Sequential_MinAdaptiveModelStep):
    def __init__(self,modelList,system=defaultSystem,stepExact=True):
        Sequential_MinAdaptiveModelStep.__init__(self,modelList,system,stepExact)
        self.modelList = modelList[:len(pnList)-1]

class Sequential_MinAdaptiveModelStepPS(Sequential_FixedStep):
    def __init__(self,modelList,system=defaultSystem,stepExact=True):
        Sequential_FixedStep.__init__(self,modelList,system,stepExact)
        self.modelList = modelList[:len(pnList)-1]

dt_system_fixed = cylinder3p.dt_fixed
systemStepControllerType = Sequential_MinAdaptiveModelStepPS

needEBQ_GLOBAL = False
needEBQ = False

systemStepExact=False

tnList = cylinder3p.tnList



info = open("TimeList.txt","w")
#archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
