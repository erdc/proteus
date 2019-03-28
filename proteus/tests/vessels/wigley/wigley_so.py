"""
The split operator module for air/water flow around a moving rigid cylinder
"""
from proteus.default_so import *
import wigley

if wigley.useOnlyVF:
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

if wigley.movingDomain:
    pnList.append(("moveMesh_p","moveMesh_n"))

if wigley.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "wigley"

#systemStepControllerType = ISO_fixed_MinAdaptiveModelStep
systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,wigley.dt_init]+[wigley.dt_init+ i*wigley.dt_out for i in range(1,wigley.nDTout+1)]

