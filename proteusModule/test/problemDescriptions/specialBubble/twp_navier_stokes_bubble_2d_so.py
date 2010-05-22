from pyadh.default_so import *
import bubble

if bubble.applyCorrection:
    pnList = [("curvature_bubble_2d_p" , 
               "curvature_bubble_2d_c0p1_n"),
              ("twp_navier_stokes_bubble_2d_p" , 
               "twp_navier_stokes_bubble_2d_c0p1c0p1_n"),
              ("ls_bubble_2d_p" , "ls_bubble_2d_c0p1_n"),
              ("vof_bubble_2d_p" , "vof_bubble_2d_c0p1_n"),
              ("redist_bubble_2d_p" , "redist_bubble_2d_c0p1_n"),
              ("ls_consrv_bubble_2d_p" , 
               "ls_consrv_bubble_2d_c0p1_n")]
elif bubble.applyRedistancing:
    pnList = [("curvature_bubble_2d_p" , 
               "curvature_bubble_2d_c0p1_n"),
              ("twp_navier_stokes_bubble_2d_p" , 
               "twp_navier_stokes_bubble_2d_c0p1c0p1_n"),
              ("ls_bubble_2d_p" , "ls_bubble_2d_c0p1_n"),
              ("redist_bubble_2d_p" , "redist_bubble_2d_c0p1_n")]
else:
    pnList = [("curvature_bubble_2d_p" , 
               "curvature_bubble_2d_c0p1_n"),
              ("twp_navier_stokes_bubble_2d_p" , 
               "twp_navier_stokes_bubble_2d_c0p1c0p1_n"),
              ("ls_bubble_2d_p" , "ls_bubble_2d_c0p1_n")]
    
name = "twp_navier_stokes_bubble_2d"
if bubble.useBackwardEuler:
    systemStepControllerType = Sequential_FixedStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
needEBQ_GLOBAL = False#True
needEBQ = False#True
useOneArchive = False
if bubble.nDTout != None:
    tnList = [i*bubble.T/(float(bubble.nDTout)-1) for i in range(bubble.nDTout)]
else:
    tnList = [0.0,1.0e-3,bubble.T]
