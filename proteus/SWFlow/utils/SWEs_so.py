"""
Split operator module for SWFlow: SWEs or Dispersive SWEs
"""
import os
from proteus.default_so import *
from proteus import Context
import sys

# (!) not the greatest for getting the file name but it works
name = str(sys.argv[0][:-3])
parun_passed = False
for i in range(len(sys.argv)):
    if 'parun' in sys.argv[i]:
        parun_passed = True
    if parun_passed is True and sys.argv[i][-3:] == '.py':
        name = sys.argv[i][:-3]
        break
    else:
        name = "SWFlow"

# ***************************** #
# ********** CONTEXT ********** #
# ***************************** #
case = __import__(name)
Context.setFromModule(case)
ct = Context.get()

# READ FROM myTpFlowProblem #
assert hasattr(
    ct, 'mySWFlowProblem'), "Create mySWFlowProblem from SWFlowProblem"
sw_model = ct.mySWFlowProblem.sw_model
outputStepping = ct.mySWFlowProblem.outputStepping

# **************************** #
# ********** pnList ********** #
# **************************** #
if sw_model == 0:  # SWEs
    pnList = [("sw_p", "sw_n")]
else: # DSWEs
    pnList = [("GN_sw_p", "GN_sw_n")]

# **************************** #
# ********** tnList ********** #
# **************************** #
tnList = [0., outputStepping['dt_init']] + [float(k) * outputStepping['final_time'] / float(
    outputStepping['nDTout']) for k in range(1, outputStepping['nDTout'] + 1)]
