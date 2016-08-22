from proteus.default_so import *
import channel

if channel.useOnlyVF:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n")]
else:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("ls_p",                "ls_n"),
              ("redist_p",            "redist_n")]
    
    
if channel.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "channel_p" 

systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,channel.dt_init]+[i*channel.dt_fixed for i in range(1,channel.nDTout+1)] 
