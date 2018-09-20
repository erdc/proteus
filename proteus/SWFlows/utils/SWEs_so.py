"""
Split operator module for SWFlow: SWEs or Dispersive SWEs
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

# READ FROM myTpFlowProblem #
assert hasattr(ct,'mySWFlowProblem'), "Create mySWFlowProblem from SWFlowProblem"
sw_model = ct.mySWFlowProblem.sw_model
outputStepping = ct.mySWFlowProblem.outputStepping

# **************************** #
# ********** pnList ********** #
# **************************** #
if sw_model==0: #SWEs
    pnList = [("sw_p", "sw_n")]
else:
    raise("Not implemented!")
    #pnList = [("dsw_p", "dsw_n")]

# **************************** #
# ********** tnList ********** #
# **************************** #
tnList=[0.,outputStepping['dt_init']]+[float(k)*outputStepping['final_time']/float(outputStepping['nDTout']) for k in range(1,outputStepping['nDTout']+1)]

