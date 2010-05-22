from pyadh.default_so import *
import wavetank3d

assert(not (wavetank3d.applyCorrection==True and wavetank3d.applyRedistancing==False))

if wavetank3d.applyCorrection and wavetank3d.applyRedistancing:
    pnList = [("twp_navier_stokes_wavetank_3d_p",
               "twp_navier_stokes_wavetank_3d_c0p1c0p1_n"),
              ("ls_wavetank_3d_p",
               "ls_wavetank_3d_c0p1_n"), 
              ("vof_wavetank_3d_p",
               "vof_wavetank_3d_c0p1_n"), 
              ("redist_wavetank_3d_p",
               "redist_wavetank_3d_c0p1_n"), 
              ("ls_consrv_wavetank_3d_p",
               "ls_consrv_wavetank_3d_c0p1_n")] 
elif wavetank.applyRedistancing:
    pnList = [("twp_navier_stokes_wavetank_3d_p",
               "twp_navier_stokes_wavetank_3d_c0p1c0p1_n"),
              ("ls_wavetank_3d_p",
               "ls_wavetank_3d_c0p1_n"), 
              ("redist_wavetank_3d_p",
               "redist_wavetank_3d_c0p1_n")]
else:
    pnList = [("twp_navier_stokes_wavetank_3d_p",
               "twp_navier_stokes_wavetank_3d_c0p1c0p1_n"),
              ("ls_wavetank_3d_p",
               "ls_wavetank_3d_c0p1_n")]

name = "wavetank"
systemStepControllerType = Sequential_MinAdaptiveModelStep
needEBQ_GLOBAL = True#False
needEBQ = True#False
useOneArchive = False
archiveFlag = ArchiveFlags.EVERY_USER_STEP
if wavetank3d.nDTout != None:
    tnList = [0.0,wavetank3d.dt_init]
    tnList.extend([wavetank3d.dt_init+i*wavetank3d.T/(float(wavetank3d.nDTout)-1) for i in range(1,wavetank3d.nDTout)])
else:
    tnList = [0.0,wavetank3d.dt_init,wavetank3d.T]

