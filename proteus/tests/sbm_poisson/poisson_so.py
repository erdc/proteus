from proteus.default_so import *
from proteus import *
import poisson2D

pnList = [("poisson_p", "poisson_n")]

systemStepControllerType = Sequential_MinAdaptiveModelStep

name = poisson2D.soname + "_ls"

needEBQ_GLOBAL = False
needEBQ = False

# archiveFlag = ArchiveFlags.EVERY_USER_STEP
#archiveFlag = ArchiveFlags.EVERY_MODEL_STEP
tnList = poisson2D.tnList
# cek hard coded steps for article snapshots
#tnList = [0.0,4.0,8.0]
useOneArchive = True

# dt_system_fixed = poisson2D.dt_fixed
# stepController = StepControl.Min_dt_cfl_controller
# systemStepControllerType  = SplitOperator.Sequential_FixedStep


systemStepExact = True

