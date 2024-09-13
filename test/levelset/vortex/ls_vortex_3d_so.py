from proteus.default_so import *
from . import vortex
from .vortex import *

if applyRedistancing:
    if applyCorrection:
        pnList = [("ls_vortex_3d_p","ls_vortex_3d_n"),
                  ("redist_vortex_3d_p","redist_vortex_3d_n"), 
                  ("vof_vortex_3d_p","vof_vortex_3d_n"), 
                  ("ls_consrv_vortex_3d_p","ls_consrv_vortex_3d_n")]
    else:
        pnList = [("ls_vortex_3d_p","ls_vortex_3d_n"),
                  ("redist_vortex_3d_p","redist_vortex_3d_n"), 
                  ("vof_vortex_3d_p","vof_vortex_3d_n")]
else:
    pnList = [("ls_vortex_3d_p","ls_vortex_3d_n"),
              ("vof_vortex_3d_p","vof_vortex_3d_n")]

if "FLCBDF" in [timeIntegration_vof,timeIntegration_ls]:
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
