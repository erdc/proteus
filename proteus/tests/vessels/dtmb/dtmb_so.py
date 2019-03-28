"""
The split operator module for air/water flow around a moving rigid cylinder
"""
from proteus.default_so import *
import dtmb

if dtmb.useOnlyVF:
    pnList = [("twp_navier_stokes_p", #0
               "twp_navier_stokes_n"),
              ("vof_p", #1              
               "vof_n")]
else:
    pnList = [("twp_navier_stokes_p" , #0
               "twp_navier_stokes_n"),
              ("vof_p" , #1
               "vof_n"),
              ("ls_p" , #2
               "ls_n"),
              ("redist_p" ,#3 
               "redist_n"),
              ("ls_consrv_p" ,#4 
               "ls_consrv_n")]

if dtmb.movingDomain:
    pnList.append(("moveMesh_p","moveMesh_n"))

if dtmb.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "dtmb"

systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,dtmb.dt_init]+[dtmb.dt_init+ i*dtmb.dt_out for i in range(1,dtmb.nDTout+1)]

