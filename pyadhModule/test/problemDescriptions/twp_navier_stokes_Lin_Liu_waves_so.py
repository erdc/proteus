from pyadh.default_so import *
import Lin_Liu_waves

name = "Lin_Liu_waves_002"
if Lin_Liu_waves.applyCorrection:
    pnList = [("curvature_Lin_Liu_waves_2d_p" , 
               "curvature_Lin_Liu_waves_2d_c0p1_n"),
              ("twp_navier_stokes_Lin_Liu_waves_2d_p" , 
               "twp_navier_stokes_Lin_Liu_waves_2d_c0p1c0p1_n"),
              ("ls_Lin_Liu_waves_2d_p" , "ls_Lin_Liu_waves_2d_c0p1_n"),
              ("vof_Lin_Liu_waves_2d_p" , "vof_Lin_Liu_waves_2d_c0p1_n"),
              ("redist_Lin_Liu_waves_2d_p" , "redist_Lin_Liu_waves_2d_c0p1_n"),
              ("ls_consrv_Lin_Liu_waves_2d_p" , 
               "ls_consrv_Lin_Liu_waves_2d_c0p1_n")]
elif Lin_Liu_waves.applyRedistancing:
    pnList = [("curvature_Lin_Liu_waves_2d_p" , 
               "curvature_Lin_Liu_waves_2d_c0p1_n"),
              ("twp_navier_stokes_Lin_Liu_waves_2d_p" , 
               "twp_navier_stokes_Lin_Liu_waves_2d_c0p1c0p1_n"),
              ("ls_Lin_Liu_waves_2d_p" , "ls_Lin_Liu_waves_2d_c0p1_n"),
              ("redist_Lin_Liu_waves_2d_p" , "redist_Lin_Liu_waves_2d_c0p1_n")]#mwf changed to c0p2
else:
    pnList = [("curvature_Lin_Liu_waves_2d_p" , 
               "curvature_Lin_Liu_waves_2d_c0p1_n"),
              ("twp_navier_stokes_Lin_Liu_waves_2d_p" , 
               "twp_navier_stokes_Lin_Liu_waves_2d_c0p1c0p1_n"),
              ("ls_Lin_Liu_waves_2d_p" , "ls_Lin_Liu_waves_2d_c0p1_n")]#mwf changed to c0p2
    
if Lin_Liu_waves.useBackwardEuler:
    #systemStepControllerType = Sequential_FixedStep
    #systemStepControllerType = Sequential_MinModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
needEBQ_GLOBAL =True# False#True
needEBQ = True#False#True
useOneArchive = False
archiveFlag = ArchiveFlags.EVERY_USER_STEP

if Lin_Liu_waves.nDTout != None:
    tnList = [0.0,Lin_Liu_waves.dt_init]
    tnList.extend([Lin_Liu_waves.dt_init+i*Lin_Liu_waves.T/(float(Lin_Liu_waves.nDTout)-1) for i in range(1,Lin_Liu_waves.nDTout)])
else:
    tnList = [0.0,Lin_Liu_waves.dt_init,Lin_Liu_waves.T]
