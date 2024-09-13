from proteus.default_so import *
try:
    from . import rotation2D
    from .rotation2D import *
except:
    import rotation2D
    from rotation2D import *

if applyRedistancing:
    if applyCorrection:
        pnList = [("ls_rotation_2d_p","ls_rotation_2d_n"),
                  ("redist_rotation_2d_p","redist_rotation_2d_n"), 
                  ("vof_rotation_2d_p","vof_rotation_2d_n"), 
                  ("ls_consrv_rotation_2d_p","ls_consrv_rotation_2d_n")]
    else:
        pnList = [("ls_rotation_2d_p","ls_rotation_2d_n"),
                  ("redist_rotation_2d_p","redist_rotation_2d_n"), 
                  ("vof_rotation_2d_p","vof_rotation_2d_n")]
elif not onlyVOF:
    pnList = [("ls_rotation_2d_p","ls_rotation_2d_n"),
              ("vof_rotation_2d_p","vof_rotation_2d_n")]
else:
    pnList = [("vof_rotation_2d_p","vof_rotation_2d_n")]

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
