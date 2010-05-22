from pyadh.default_so import *
import obstacleInFlume3d

if obstacleInFlume3d.applyCorrection:
    pnList = [("twp_navier_stokes_obstacleInFlume_3d_p" , 
               "twp_navier_stokes_obstacleInFlume_3d_n"),
              ("ls_obstacleInFlume_3d_p" , "ls_obstacleInFlume_3d_n"),
              ("vof_obstacleInFlume_3d_p" , "vof_obstacleInFlume_3d_n"),
              ("redist_obstacleInFlume_3d_p" , "redist_obstacleInFlume_3d_n"),
              ("ls_consrv_obstacleInFlume_3d_p" , 
               "ls_consrv_obstacleInFlume_3d_n")]
elif obstacleInFlume3d.applyRedistancing:
    pnList = [("twp_navier_stokes_obstacleInFlume_3d_p" , 
               "twp_navier_stokes_obstacleInFlume_3d_n"),
              ("ls_obstacleInFlume_3d_p" , "ls_obstacleInFlume_3d_n"),
              ("redist_obstacleInFlume_3d_p" , "redist_obstacleInFlume_3d_n")]
else:
    pnList = [("twp_navier_stokes_obstacleInFlume_3d_p" , 
               "twp_navier_stokes_obstacleInFlume_3d_n"),
              ("ls_obstacleInFlume_3d_p" , "ls_obstacleInFlume_3d_n")]
    
name = "twp_navier_stokes_obstacleInFlume_3d"

if obstacleInFlume3d.useBackwardEuler:
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
needEBQ_GLOBAL = False#True
needEBQ = False#True
useOneArchive = False

tnList = [0.0,obstacleInFlume3d.dt_init]+[obstacleInFlume3d.dt_init+i*(obstacleInFlume3d.T-obstacleInFlume3d.dt_init)/float(obstacleInFlume3d.nDTout-1) for i in range(1,obstacleInFlume3d.nDTout)]
