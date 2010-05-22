from pyadh.default_so import *
import beach_erosion_board_waves_3d

name = "beach_erosion_board_waves_001"
if beach_erosion_board_waves_3d.applyCorrection:
    pnList = [("twp_navier_stokes_beach_erosion_board_waves_3d_p" , 
               "twp_navier_stokes_beach_erosion_board_waves_3d_c0p1c0p1_n"),
              ("ls_beach_erosion_board_waves_3d_p" , "ls_beach_erosion_board_waves_3d_c0p1_n"),
              ("vof_beach_erosion_board_waves_3d_p" , "vof_beach_erosion_board_waves_3d_c0p1_n"),
              ("redist_beach_erosion_board_waves_3d_p" , "redist_beach_erosion_board_waves_3d_c0p1_n"),
              ("ls_consrv_beach_erosion_board_waves_3d_p" , 
               "ls_consrv_beach_erosion_board_waves_3d_c0p1_n")]
elif beach_erosion_board_waves_3d.applyRedistancing:
    pnList = [("twp_navier_stokes_beach_erosion_board_waves_3d_p" , 
               "twp_navier_stokes_beach_erosion_board_waves_3d_c0p1c0p1_n"),
              ("ls_beach_erosion_board_waves_3d_p" , "ls_beach_erosion_board_waves_3d_c0p1_n"),
              ("redist_beach_erosion_board_waves_3d_p" , "redist_beach_erosion_board_waves_3d_c0p1_n")]
else:
    pnList = [("twp_navier_stokes_beach_erosion_board_waves_3d_p" , 
               "twp_navier_stokes_beach_erosion_board_waves_3d_c0p1c0p1_n"),
              ("ls_beach_erosion_board_waves_3d_p" , "ls_beach_erosion_board_waves_3d_c0p1_n")]
    
if beach_erosion_board_waves_3d.useBackwardEuler:
    #systemStepControllerType = Sequential_FixedStep
    #systemStepControllerType = Sequential_MinModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
needEBQ_GLOBAL =True
needEBQ = True
useOneArchive = False
archiveFlag = ArchiveFlags.EVERY_USER_STEP

if beach_erosion_board_waves_3d.nDTout != None:
    tnList = [0.0,beach_erosion_board_waves_3d.dt_init]
    tnList.extend([beach_erosion_board_waves_3d.dt_init+i*beach_erosion_board_waves_3d.T/(float(beach_erosion_board_waves_3d.nDTout)-1) for i in range(1,beach_erosion_board_waves_3d.nDTout)])
else:
    tnList = [0.0,beach_erosion_board_waves_3d.dt_init,beach_erosion_board_waves_3d.T]
