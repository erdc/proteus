"""
Split operator module for multiphase: CLSVOF with RANS2P
"""
import os
from proteus.default_so import *
from . import multiphase

if multiphase.useCLSVOF:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("clsvof_p",               "clsvof_n")]
else:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n"),
              ("ls_p",                "ls_n"),
              ("redist_p",            "redist_n"),
              ("ls_consrv_p",         "ls_consrv_n")]
    
name = "multiphase" 

if multiphase.timeDiscretization == 'flcbdf':
    systemStepControllerType = Sequential_MinFLCBDFModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,multiphase.dt_init]+[i*multiphase.dt_fixed for i in range(1,multiphase.nDTout+1)]

info = open("TimeList.txt","w")

