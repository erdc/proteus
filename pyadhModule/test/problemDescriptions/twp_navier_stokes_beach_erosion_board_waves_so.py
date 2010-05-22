from pyadh.default_so import *
import beach_erosion_board_waves

name = "beach_erosion_board_waves_002"
if beach_erosion_board_waves.applyCorrection:
    pnList = [("twp_navier_stokes_beach_erosion_board_waves_2d_p" , 
               "twp_navier_stokes_beach_erosion_board_waves_2d_c0p1c0p1_n"),
              ("ls_beach_erosion_board_waves_2d_p" , "ls_beach_erosion_board_waves_2d_c0p1_n"),
              ("vof_beach_erosion_board_waves_2d_p" , "vof_beach_erosion_board_waves_2d_c0p1_n"),
              ("redist_beach_erosion_board_waves_2d_p" , "redist_beach_erosion_board_waves_2d_c0p1_n"),
              ("ls_consrv_beach_erosion_board_waves_2d_p" , 
               "ls_consrv_beach_erosion_board_waves_2d_c0p1_n")]
elif beach_erosion_board_waves.applyRedistancing:
    pnList = [("twp_navier_stokes_beach_erosion_board_waves_2d_p" , 
               "twp_navier_stokes_beach_erosion_board_waves_2d_c0p1c0p1_n"),
              ("ls_beach_erosion_board_waves_2d_p" , "ls_beach_erosion_board_waves_2d_c0p1_n"),
              ("redist_beach_erosion_board_waves_2d_p" , "redist_beach_erosion_board_waves_2d_c0p1_n")]#mwf changed to c0p2
else:
    pnList = [("twp_navier_stokes_beach_erosion_board_waves_2d_p" , 
               "twp_navier_stokes_beach_erosion_board_waves_2d_c0p1c0p1_n"),
              ("ls_beach_erosion_board_waves_2d_p" , "ls_beach_erosion_board_waves_2d_c0p1_n")]#mwf changed to c0p2
    
if beach_erosion_board_waves.useBackwardEuler:
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinFLCBDFModelStep
needEBQ_GLOBAL =True# False#True
needEBQ = True#False#True
useOneArchive = False
archiveFlag = ArchiveFlags.EVERY_USER_STEP

if beach_erosion_board_waves.nDTout != None:
    tnList = [0.0,beach_erosion_board_waves.dt_init]
    tnList.extend([beach_erosion_board_waves.dt_init+i*beach_erosion_board_waves.T/(float(beach_erosion_board_waves.nDTout)-1) for i in range(1,beach_erosion_board_waves.nDTout)])
else:
    tnList = [0.0,beach_erosion_board_waves.dt_init,beach_erosion_board_waves.T]
