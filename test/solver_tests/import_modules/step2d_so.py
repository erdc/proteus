from proteus.default_so import *
from proteus import defaults
defaults.reset_default_so()
import proteus
try:
    from . import step2d
except:
    import step2d

from proteus import Context
import os

Context.setFromModule(step2d)
ct = Context.get()

pnList = [("twp_navier_stokes_step2d_p", "twp_navier_stokes_step2d_n")]

name = "twp_navier_stokes_step2d"

systemStepControllerType = proteus.SplitOperator.Sequential_tnList
systemStepExact = True

tnList = [0.0, 1.0]
