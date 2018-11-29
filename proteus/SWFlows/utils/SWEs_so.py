"""
Split operator module for SWFlow: SWEs or Dispersive SWEs
"""
from __future__ import absolute_import
from builtins import range
import os
from proteus.default_so import *
from proteus import Context
import sys

# (!) not the greatest for getting the file name but it works
#name = str(sys.argv[0][:-3])
name = "SWFlow"
for i in range(len(sys.argv)):
    if '-f' in sys.argv[i]:
        assert sys.argv[i+1][-3:], "fileName must end with .py"
        name = sys.argv[i+1][:-3]
        break

# ***************************** #
# ********** CONTEXT ********** #
# ***************************** #
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
    pnList = [("dsw_p", "dsw_n")]
    #raise("Not implemented!")

# **************************** #
# ********** tnList ********** #
# **************************** #
tnList=[0.,outputStepping['dt_init']]+[float(k)*outputStepping['final_time']/float(outputStepping['nDTout']) for k in range(1,outputStepping['nDTout']+1)]
