from pyadh.default_so import *
import damBreak

if damBreak.applyCorrection:
    pnList = [("curvature_damBreak_2d_p" , 
               "curvature_damBreak_2d_n"),
              ("twp_navier_stokes_damBreak_2d_p" , 
               "twp_navier_stokes_damBreak_2d_n"),
              ("ls_damBreak_2d_p" , "ls_damBreak_2d_n"),
              ("vof_damBreak_2d_p" , "vof_damBreak_2d_n"),
              ("redist_damBreak_2d_p" , "redist_damBreak_2d_n"),
              ("ls_consrv_damBreak_2d_p" , 
               "ls_consrv_damBreak_2d_n")]
elif damBreak.applyRedistancing:
    pnList = [("curvature_damBreak_2d_p" , 
               "curvature_damBreak_2d_n"),
              ("twp_navier_stokes_damBreak_2d_p" , 
               "twp_navier_stokes_damBreak_2d_n"),
              ("ls_damBreak_2d_p" , "ls_damBreak_2d_n"),
              ("redist_damBreak_2d_p" , "redist_damBreak_2d_n")]
else:
    pnList = [("curvature_damBreak_2d_p" , 
               "curvature_damBreak_2d_n"),
              ("twp_navier_stokes_damBreak_2d_p" , 
               "twp_navier_stokes_damBreak_2d_n"),
              ("ls_damBreak_2d_p" , "ls_damBreak_2d_n")]
    
name = "twp_navier_stokes_damBreak_2d"
if damBreak.useBackwardEuler:
    systemStepControllerType = Sequential_FixedStep
    systemStepControllerType = Sequential_MinModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
needEBQ_GLOBAL = True
needEBQ = True
useOneArchive = False
archiveFlag = ArchiveFlags.EVERY_USER_STEP
if damBreak.nDTout != None:
    tnList = [0.0,damBreak.dt_init]+[damBreak.dt_init+i*(damBreak.T-damBreak.dt_init)/(float(damBreak.nDTout)-1) for i in range(1,damBreak.nDTout)]
else:
    tnList = [0.0,damBreak.dt_init,damBreak.T]
