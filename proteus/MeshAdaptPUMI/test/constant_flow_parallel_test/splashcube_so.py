from proteus.default_so import *
import splashcube
from proteus.MeshAdaptPUMI import MeshAdaptPUMI

if splashcube.useOnlyVF:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n")]
else:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n"),
              ("ls_p",                "ls_n"),
              ("redist_p",            "redist_n"),
              ("ls_consrv_p",         "ls_consrv_n")]


if splashcube.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "splashcube_p"

systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,splashcube.dt_init,2.0*splashcube.dt_init]+[i*splashcube.dt_fixed for i in range(1,splashcube.nDTout+1)]
archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
