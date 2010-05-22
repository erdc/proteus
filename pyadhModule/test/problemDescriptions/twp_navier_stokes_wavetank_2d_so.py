from pyadh.default_so import *
import wavetank

assert(not (wavetank.applyCorrection==True and wavetank.applyRedistancing==False))

if wavetank.applyCorrection and wavetank.applyRedistancing:
    pnList = [("twp_navier_stokes_wavetank_2d_p",
               "twp_navier_stokes_wavetank_2d_c0p1c0p1_n"),
              ("ls_wavetank_2d_p",
               "ls_wavetank_2d_c0p1_n"), 
              ("vof_wavetank_2d_p",
               "vof_wavetank_2d_c0p1_n"), 
              ("redist_wavetank_2d_p",
               "redist_wavetank_2d_c0p1_n"), 
              ("ls_consrv_wavetank_2d_p",
               "ls_consrv_wavetank_2d_c0p1_n")] 
elif wavetank.applyRedistancing:
    pnList = [("twp_navier_stokes_wavetank_2d_p",
               "twp_navier_stokes_wavetank_2d_c0p1c0p1_n"),
              ("ls_wavetank_2d_p",
               "ls_wavetank_2d_c0p1_n"), 
              ("redist_wavetank_2d_p",
               "redist_wavetank_2d_c0p1_n")]
else:
    pnList = [("twp_navier_stokes_wavetank_2d_p",
               "twp_navier_stokes_wavetank_2d_c0p1c0p1_n"),
              ("ls_wavetank_2d_p",
               "ls_wavetank_2d_c0p1_n")]

name = "wavetank"
systemStepControllerType = Sequential_MinAdaptiveModelStep
needEBQ_GLOBAL = True#False
needEBQ = True#False
useOneArchive = False
archiveFlag = ArchiveFlags.EVERY_USER_STEP
tnList = [0.0,wavetank.dt_init,wavetank.T]

