from proteus.default_so import *
from proteus import defaults
defaults.reset_default_so()
import proteus
try:
    from . import cavity2d
except:
    import cavity2d

from proteus import Context
import os

Context.setFromModule(cavity2d)
ct = Context.get()

pnList = [("twp_navier_stokes_cylinder_2d_p", "twp_navier_stokes_cylinder_2d_n")]

name = "twp_navier_stokes_cavity_2d"

systemStepControllerType = proteus.SplitOperator.Sequential_tnList
systemStepExact = True

tnList = [0.0,100000.0]
