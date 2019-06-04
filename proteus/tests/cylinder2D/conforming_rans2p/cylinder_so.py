from proteus import *
import proteus.default_so
reload(proteus.default_so)
from proteus.default_so import *
try:
    from . import cylinder2d as cylinder
except:
    import cylinder2d as cylinder

from proteus.SplitOperator import Sequential_FixedStep_Simple, defaultSystem

pnList = [("twp_navier_stokes_cylinder_2d_p", "twp_navier_stokes_cylinder_2d_n")]
name = "cylinder"
soname=name

class Sequential_MinAdaptiveModelStepPS(Sequential_FixedStep):
    def __init__(self,modelList,system=defaultSystem,stepExact=True):
        Sequential_FixedStep.__init__(self,modelList,system,stepExact)
        self.modelList = modelList

dt_system_fixed = cylinder.dt_fixed
#systemStepControllerType = Sequential_MinAdaptiveModelStepPS

systemStepControllerType = Sequential_FixedStep #Sequential_FixedStep_Simple # uses time steps in so.tnList
# dt_system_fixed = 0.01; 
systemStepExact=False;
# 
# 
# needEBQ_GLOBAL = False
# needEBQ = False
# 
# tnList = [0.0,cylinder.dt_init]+[i*cylinder.dt_fixed for i in range(1,cylinder.nDTout+1)]
tnList = cylinder.tnList
#archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
