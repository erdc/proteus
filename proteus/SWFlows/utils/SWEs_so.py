"""
Split operator module for multiphase: CLSVOF with RANS2P
"""
from __future__ import absolute_import
from builtins import range
import os
from proteus.default_so import *
from proteus import Context

# ***************************** #
# ********** CONTEXT ********** #
# ***************************** #
name = "SWFlow"
case = __import__(name)
Context.setFromModule(case)
ct = Context.get()

# **************************** #
# ********** pnList ********** #
# **************************** #
pnList = [("sw_p", "sw_n")]


tnList=[0.,1E-6]+[float(n)*ct.T/float(ct.nDTout) for n in range(1,ct.nDTout+1)]
