from __future__ import absolute_import
from proteus.default_so import *
import proteus
from . import step2d
reload(step2d)
from proteus import Context
import os

Context.setFromModule(step2d)
ct = Context.get()

pnList = [("twp_navier_stokes_step2d_p", "twp_navier_stokes_step2d_n")]

name = "twp_navier_stokes_cylinder_2d"

systemStepControllerType = proteus.SplitOperator.Sequential_tnList
systemStepExact = True

tnList = [0.0, 1000000000000000000.0] #, 5, 7, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 1000] #, 100]
#tnList = [0.0,cavity2d.dt_init]+[i*cavity2d.dt_fixed for i in range(1,cavity2d.nDTout+1)] 
