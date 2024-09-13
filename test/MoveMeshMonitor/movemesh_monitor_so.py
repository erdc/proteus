"""
Split operator module for two-phase flow
"""

import os
from proteus.default_so import *
from proteus import Context
try:
    from . import movemesh_monitor
except:
    import movemesh_monitor
# Create context from main module
name_so = os.path.basename(__file__)
if '_so.py' in name_so[-6:]:
    name = name_so[:-6]
elif '_so.pyc' in name_so[-7:]:
    name = name_so[:-7]
else:
    raise NameError('Split operator module must end with "_so.py"')

Context.setFromModule(movemesh_monitor)
ct = Context.get()

# List of p/n files
pnList = []

# moving mesh
if ct.movingDomain:
    pnList += [("moveMesh_p", "moveMesh_n")]

tnList = [0., ct.dt_init, ct.T]
dt_system_fixed = 0.001
systemStepControllerType = Sequential_FixedStep
needEBQ_GLOBAL = False
needEBQ = False

from proteus.Archiver import ArchiveFlags
archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
