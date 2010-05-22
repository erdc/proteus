from pyadh.default_so import *
import sloshbox3d

if sloshbox3d.applyCorrection:
    pnList = [("twp_navier_stokes_sloshbox_3d_old_p" , 
               "twp_navier_stokes_sloshbox_3d_c0p1c0p1_old_n"),
              ("ls_sloshbox_3d_old_p" , "ls_sloshbox_3d_c0p1_old_n"),
              ("vof_sloshbox_3d_old_p" , "vof_sloshbox_3d_c0p1_old_n"),
              ("redist_sloshbox_3d_old_p" , "redist_sloshbox_3d_c0p1_old_n"),
              ("ls_consrv_sloshbox_3d_old_p" , 
               "ls_consrv_sloshbox_3d_c0p1_old_n")]
elif sloshbox3d.applyRedistancing:
    pnList = [("twp_navier_stokes_sloshbox_3d_old_p" , 
               "twp_navier_stokes_sloshbox_3d_c0p1c0p1_old_n"),
              ("ls_sloshbox_3d_old_p" , "ls_sloshbox_3d_c0p1_old_n"),
              ("redist_sloshbox_3d_old_p" , "redist_sloshbox_3d_c0p1_old_n")]
else:
    pnList = [("twp_navier_stokes_sloshbox_3d_old_p" , 
               "twp_navier_stokes_sloshbox_3d_c0p1c0p1_old_n"),
              ("ls_sloshbox_3d_old_p" , "ls_sloshbox_3d_c0p1_old_n")]
    
name = "twp_navier_stokes_sloshbox_3d"

if sloshbox3d.useBackwardEuler:
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
needEBQ_GLOBAL = True
needEBQ = True
useOneArchive = False

tnList = [0.0,sloshbox3d.dt_init]+[sloshbox3d.dt_init+i*(sloshbox3d.T-sloshbox3d.dt_init)/float(sloshbox3d.nDTout-1) for i in range(1,sloshbox3d.nDTout)]
