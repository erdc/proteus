from proteus.default_so import *
import dambreak
from proteus import MeshAdaptPUMI

if dambreak.useOnlyVF:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n")]
else:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n"),
              ("ls_p",                "ls_n"),
              ("redist_p",            "redist_n"),
              ("ls_consrv_p",         "ls_consrv_n")]


if dambreak.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "dambreak_p"

#systemStepControllerType = Sequential_MinAdaptiveModelStep
systemStepControllerType = Sequential_FixedStep_Simple

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,dambreak.dt_init]+[i*dambreak.dt_fixed for i in range(1,dambreak.nDTout+1)]
archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
#archiveFlag = ArchiveFlags.EVERY_USER_STEP
