from pyadh.default_so import *
import embankmentRunup

name = "cobras_embankment_002"
if embankmentRunup.applyCorrection:
    pnList = [("curvature_cobras_embankment_2d_p" , 
               "curvature_cobras_embankment_2d_c0p1_n"),
              ("twp_navier_stokes_cobras_embankment_2d_p" , 
               "twp_navier_stokes_cobras_embankment_2d_c0p1c0p1_n"),
              ("ls_cobras_embankment_2d_p" , "ls_cobras_embankment_2d_c0p1_n"),
              ("vof_cobras_embankment_2d_p" , "vof_cobras_embankment_2d_c0p1_n"),
              ("redist_cobras_embankment_2d_p" , "redist_cobras_embankment_2d_c0p1_n"),
              ("ls_consrv_cobras_embankment_2d_p" , 
               "ls_consrv_cobras_embankment_2d_c0p1_n")]
elif embankmentRunup.applyRedistancing:
    pnList = [("curvature_cobras_embankment_2d_p" , 
               "curvature_cobras_embankment_2d_c0p1_n"),
              ("twp_navier_stokes_cobras_embankment_2d_p" , 
               "twp_navier_stokes_cobras_embankment_2d_c0p1c0p1_n"),
              ("ls_cobras_embankment_2d_p" , "ls_cobras_embankment_2d_c0p1_n"),
              ("redist_cobras_embankment_2d_p" , "redist_cobras_embankment_2d_c0p1_n")]#mwf changed to c0p2
else:
    pnList = [("curvature_cobras_embankment_2d_p" , 
               "curvature_cobras_embankment_2d_c0p1_n"),
              ("twp_navier_stokes_cobras_embankment_2d_p" , 
               "twp_navier_stokes_cobras_embankment_2d_c0p1c0p1_n"),
              ("ls_cobras_embankment_2d_p" , "ls_cobras_embankment_2d_c0p1_n")]#mwf changed to c0p2
    
if embankmentRunup.useBackwardEuler:
    #systemStepControllerType = Sequential_FixedStep
    #systemStepControllerType = Sequential_MinModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
needEBQ_GLOBAL =True# False#True
needEBQ = True#False#True
useOneArchive = False
archiveFlag = ArchiveFlags.EVERY_USER_STEP

if embankmentRunup.nDTout != None:
    tnList = [0.0,embankmentRunup.dt_init]
    tnList.extend([embankmentRunup.dt_init+i*embankmentRunup.T/(float(embankmentRunup.nDTout)-1) for i in range(1,embankmentRunup.nDTout)])
else:
    tnList = [0.0,embankmentRunup.dt_init,embankmentRunup.T]
