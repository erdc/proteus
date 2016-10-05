from proteus.default_so import *
import rotation2D
from rotation2D import *

pnList = [("ls_rotation_2d_p","ls_rotation_2d_n")]

if "FLCBDF" in [timeIntegration_ls]:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep

systemStepExact = True

name=soname

needEBQ_GLOBAL  = False
needEBQ = False

archiveFlag = ArchiveFlags.EVERY_USER_STEP
#archiveFlag = ArchiveFlags.EVERY_MODEL_STEP
DT = T/float(nDTout)
tnList = [i*DT for i  in range(nDTout+1)]
#cek hard coded steps for article snapshots
#tnList = [0.0,4.0,8.0]
useOneArchive = True
