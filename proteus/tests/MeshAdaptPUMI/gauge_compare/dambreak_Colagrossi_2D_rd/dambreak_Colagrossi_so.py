"""
Split operator module for two-phase flow
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


# List of p/n files
pnList = []

if ct.useOnlyVF:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p", "vof_n")]
else:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p", "vof_n"),
              ("ls_p", "ls_n"),
              ("redist_p", "redist_n"),
              ("ls_consrv_p", "ls_consrv_n")
             ]

if ct.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "dambreak_Colagrossi_p"

if ct.timeDiscretization == 'flcbdf':
    systemStepControllerType = Sequential_MinFLCBDFModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep

#systemStepControllerType = Sequential_FixedStep #Sequential_FixedStep_Simple # uses time steps in so.tnList
if(ct.fixedTimeStep):
    systemStepControllerType = Sequential_FixedStep_Simple # uses time steps in so.tnList
#dt_system_fixed = 0.01; 
systemStepExact=False;


needEBQ_GLOBAL = False
needEBQ = False

tnList=[0.0,ct.dt_init]+[ct.dt_init+ i*ct.dt_fixed for i in range(1,ct.nDTout+1)]