"""
Split operator module for multiphase: CLSVOF with RANS2P
"""
from __future__ import absolute_import
from builtins import range
import os
from proteus.default_so import *
from proteus import Context

# ASSUME CLSVOF #
pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
          ("clsvof_p",               "clsvof_n")]

name = "TwoPhaseFlow"
case = __import__(name)
Context.setFromModule(case)
ct = Context.get()

systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,ct.dt_init]+[i*ct.dt_output for i in range(1,ct.nDTout+1)]
info = open("TimeList.txt","w")
