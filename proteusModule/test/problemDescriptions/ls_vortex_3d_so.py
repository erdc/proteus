from pyadh.default_so import *
import vortex
from vortex import *

if applyRedistancing:
    if applyCorrection:
        pnList = [("ls_vortex_3d_p","ls_vortex_3d_c0p1_n"),
                  ("redist_vortex_3d_p","redist_vortex_3d_c0p1_n"), 
                  ("vof_vortex_3d_p","vof_vortex_3d_c0p1_n"), 
                  ("ls_consrv_vortex_3d_p","ls_consrv_vortex_3d_c0p1_n")]
    else:
        pnList = [("ls_vortex_3d_p","ls_vortex_3d_c0p1_n"),
                  ("redist_vortex_3d_p","redist_vortex_3d_c0p1_n"), 
                  ("vof_vortex_3d_p","vof_vortex_3d_c0p1_n")]
else:
    pnList = [("ls_vortex_3d_p","ls_vortex_3d_c0p1_n"),
              ("vof_vortex_3d_p","vof_vortex_3d_c0p1_n")]
if "FLCBDF" in [timeIntegration_vof,timeIntegration_ls]:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
    #systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep
#pnList = [("ls_vortex_2d_p","ls_vortex_2d_c0p1_n")]
#name="vortex_2d_c0p1nn161"
#pnList = [("ls_vortex_2d_p","ls_vortex_2d_dgp2_n")]#also dgp1,dgpk
name=soname
if tryOpt:
    needEBQ_GLOBAL  = False
    needEBQ = False
else:
    needEBQ_GLOBAL  = True
    needEBQ = True
if tryOpt:
    nDTout = 21
    archiveFlag = ArchiveFlags.EVERY_USER_STEP
    
else:
    nDTout = 1
DT = T/float(nDTout)
tnList = [i*DT for i  in range(nDTout+1)]
#cek hard coded steps for article snapshots
tnList = [0.0,4.0,8.0]
useOneArchive = True
