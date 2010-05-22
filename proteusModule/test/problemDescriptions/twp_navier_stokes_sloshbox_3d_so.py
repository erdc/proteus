from pyadh.default_so import *
import sloshbox3d
from sloshbox3d import *

if sloshbox3d.applyCorrection:
    pnList = [("twp_navier_stokes_sloshbox_3d_p" , 
               "twp_navier_stokes_sloshbox_3d_c0p1c0p1_n"),
              ("ls_sloshbox_3d_p" , "ls_sloshbox_3d_c0p1_n"),
              ("vof_sloshbox_3d_p" , "vof_sloshbox_3d_c0p1_n"),
              ("redist_sloshbox_3d_p" , "redist_sloshbox_3d_c0p1_n"),
              ("ls_consrv_sloshbox_3d_p" , 
               "ls_consrv_sloshbox_3d_c0p1_n")]
elif sloshbox3d.applyRedistancing:
    pnList = [("twp_navier_stokes_sloshbox_3d_p" , 
               "twp_navier_stokes_sloshbox_3d_c0p1c0p1_n"),
              ("ls_sloshbox_3d_p" , "ls_sloshbox_3d_c0p1_n"),
              ("redist_sloshbox_3d_p" , "redist_sloshbox_3d_c0p1_n")]
else:
    pnList = [("twp_navier_stokes_sloshbox_3d_p" , 
               "twp_navier_stokes_sloshbox_3d_c0p1c0p1_n"),
              ("ls_sloshbox_3d_p" , "ls_sloshbox_3d_c0p1_n")]
    
name = "twp_navier_stokes_sloshbox_3d"

if sloshbox3d.useBackwardEuler:
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
needEBQ_GLOBAL = False#True
needEBQ = False#True
useOneArchive = True

tnList = [0.0,sloshbox3d.dt_init]+[sloshbox3d.dt_init+i*(sloshbox3d.T-sloshbox3d.dt_init)/float(sloshbox3d.nDTout-1) for i in range(1,sloshbox3d.nDTout)]
tnList = [0.0,sloshbox3d.dt_init,1.0,2.0,10.0,20.0]
