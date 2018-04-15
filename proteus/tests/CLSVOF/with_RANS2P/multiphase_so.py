"""
Split operator module for multiphase: CLSVOF with RANS2P
"""
import os
from proteus.default_so import *
from proteus import Context

# Create context from main module
name_so = os.path.basename(__file__)
if '_so.py' in name_so[-6:]:
    name = name_so[:-6]
elif '_so.pyc' in name_so[-7:]:
    name = name_so[:-7]
else:
    raise NameError('Split operator module must end with "_so.py"')

try:
    case = __import__(name)
    Context.setFromModule(case)
    ct = Context.get()
except ImportError:
    raise ImportError(str(name) + '.py not found')

if ct.useCLSVOF:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("clsvof_p",               "clsvof_n")]
else:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n"),
              ("ls_p",                "ls_n"),
              ("redist_p",            "redist_n"),
              ("ls_consrv_p",         "ls_consrv_n")]
    
name = "multiphase" 

if ct.timeDiscretization == 'flcbdf':
    systemStepControllerType = Sequential_MinFLCBDFModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,ct.dt_init]+[i*ct.dt_fixed for i in range(1,ct.nDTout+1)]

info = open("TimeList.txt","w")


