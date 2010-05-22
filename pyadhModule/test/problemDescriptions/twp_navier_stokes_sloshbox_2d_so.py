from pyadh.default_so import *
import sloshbox
from sloshbox import *
if sloshbox.applyCorrection:
    pnList = [("twp_navier_stokes_sloshbox_2d_p" , 
               "twp_navier_stokes_sloshbox_2d_c0p1c0p1_n"),
              ("ls_sloshbox_2d_p" , "ls_sloshbox_2d_c0p1_n"),
              ("vof_sloshbox_2d_p" , "vof_sloshbox_2d_c0p1_n"),
              ("redist_sloshbox_2d_p" , "redist_sloshbox_2d_c0p1_n"),
              ("ls_consrv_sloshbox_2d_p" , 
               "ls_consrv_sloshbox_2d_c0p1_n")]
elif sloshbox.applyRedistancing:
    pnList = [("twp_navier_stokes_sloshbox_2d_p" , 
               "twp_navier_stokes_sloshbox_2d_c0p1c0p1_n"),
              ("ls_sloshbox_2d_p" , "ls_sloshbox_2d_c0p1_n"),
              ("redist_sloshbox_2d_p" , "redist_sloshbox_2d_c0p1_n")]
else:
    pnList = [("twp_navier_stokes_sloshbox_2d_p" , 
               "twp_navier_stokes_sloshbox_2d_c0p1c0p1_n"),
              ("ls_sloshbox_2d_p" , "ls_sloshbox_2d_c0p1_n")]
    
name = soname#"twp_navier_stokes_sloshbox_2d"
if (sloshbox.useBackwardEuler and sloshbox.useBackwardEuler_ls and sloshbox.useBackwardEuler_vof):
    systemStepControllerType = Sequential_FixedStep
    systemStepControllerType = Sequential_MinModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
needEBQ_GLOBAL = True
needEBQ = True
useOneArchive = True
archiveFlag = ArchiveFlags.EVERY_USER_STEP
if sloshbox.nDTout != None:
    tnList = [0.0,sloshbox.dt_init]+[sloshbox.dt_init+i*(sloshbox.T-sloshbox.dt_init)/(float(sloshbox.nDTout)-1) for i in range(1,sloshbox.nDTout)]
else:
    tnList = [0.0,sloshbox.dt_init,sloshbox.T]
