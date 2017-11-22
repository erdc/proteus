from proteus.default_so import *
import burgers2D
from burgers2D import *

pnList = [("ls_burgers_2d_p", "ls_burgers_2d_n")]

if "FLCBDF" == ct.timeIntegration_ls:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep

systemStepExact = True

name = soname + "_ls"

needEBQ_GLOBAL = False
needEBQ = False

archiveFlag = ArchiveFlags.EVERY_USER_STEP
#archiveFlag = ArchiveFlags.EVERY_MODEL_STEP
DT = T / float(nDTout)
tnList = [i * DT for i in range(nDTout + 1)]
# cek hard coded steps for article snapshots
#tnList = [0.0,4.0,8.0]
useOneArchive = True
