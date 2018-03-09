from proteus.default_so import *
import proteus
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
