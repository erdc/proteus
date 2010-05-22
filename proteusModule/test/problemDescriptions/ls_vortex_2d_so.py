from pyadh.default_so import *
import vortex
from vortex import *

if applyRedistancing:
    if applyCorrection:
        pnList = [("ls_vortex_2d_p","ls_vortex_2d_c0p1_n"),
                  ("redist_vortex_2d_p","redist_vortex_2d_c0p1_n"), 
                  ("vof_vortex_2d_p","vof_vortex_2d_c0p1_n"), 
                  ("ls_consrv_vortex_2d_p","ls_consrv_vortex_2d_c0p1_n")]
    else:
        pnList = [("ls_vortex_2d_p","ls_vortex_2d_c0p1_n"),
                  ("redist_vortex_2d_p","redist_vortex_2d_c0p1_n"), 
                  ("vof_vortex_2d_p","vof_vortex_2d_c0p1_n")]
else:
    pnList = [("ls_vortex_2d_p","ls_vortex_2d_c0p1_n"),
              ("vof_vortex_2d_p","vof_vortex_2d_c0p1_n")]

if "FLCBDF" in [timeIntegration_vof,timeIntegration_ls]:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep
name=soname
needEBQ_GLOBAL  = False
needEBQ = False
archiveFlag = ArchiveFlags.EVERY_USER_STEP
DT = T/float(nDTout)
tnList = [i*DT for i  in range(nDTout+1)]
#cek hard coded steps for article snapshots
if abs(T-8.0) <= 1.0e-6:
    tnList = [0.0,1.0,2.0,4.0,6.0,7.0,8.0]
useOneArchive = True
