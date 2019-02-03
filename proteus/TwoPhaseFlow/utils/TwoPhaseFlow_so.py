"""
Split operator module for multiphase: CLSVOF with RANS2P
"""
from __future__ import absolute_import
from builtins import range
import os
from proteus.default_so import *
from proteus import Context
import sys

# ***************************** #
# ********** CONTEXT ********** #
# ***************************** #
name = "TwoPhaseFlow"
pathMyTpFlowProblem = None

for i in range(len(sys.argv)):
    if '-f' == sys.argv[i]:
        assert sys.argv[i+1][-3:], "fileName must end with .py"
        name = sys.argv[i+1][:-3]
    if '--path' == sys.argv[i]:
        pathMyTpFlowProblem = sys.argv[i+1]

# Load module
if pathMyTpFlowProblem is not None:
    sys.path.append(pathMyTpFlowProblem)
case = __import__(name)

# Create context
Context.setFromModule(case)
ct = Context.get()


# READ FROM myTpFlowProblem #
assert hasattr(ct,'myTpFlowProblem'), "Create myTpFlowProblem from TwoPhaseFlowProblem"
ns_model = ct.myTpFlowProblem.ns_model
outputStepping = ct.myTpFlowProblem.outputStepping

# **************************** #
# ********** pnList ********** #
# **************************** #
# ASSUME CLSVOF #
if ns_model==0: #rans2p
    pnList = [("rans2p_p", "rans2p_n"),
              ("clsvof_p", "clsvof_n")]
else: #rans3p
    pnList = [("clsvof_p", "clsvof_n"),#0
              ("rans3p_p", "rans3p_n"),#1
              ("pressureincrement_p", "pressureincrement_n"),#2
              ("pressure_p", "pressure_n"),#3
              ("pressureInitial_p", "pressureInitial_n")]#4

# ****************************************** #
# ********** TIME STEP CONTROLLER ********** #
# ****************************************** #
systemStepControllerType = Sequential_MinAdaptiveModelStep
if ns_model==1: #rans3p
    PINIT_model=4
    modelSpinUpList = [PINIT_model]
    class Sequential_MinAdaptiveModelStepPS(Sequential_MinAdaptiveModelStep):
        def __init__(self,modelList,system=defaultSystem,stepExact=True):
            Sequential_MinAdaptiveModelStep.__init__(self,modelList,system,stepExact)
            self.modelList = modelList[:len(pnList)-1]
    #
    systemStepControllerType = Sequential_MinAdaptiveModelStepPS
#
needEBQ_GLOBAL = False
needEBQ = False

# **************************** #
# ********** tnList ********** #
# **************************** #
tnList=[0.,outputStepping['dt_init']]+[float(k)*outputStepping['final_time']/float(outputStepping['nDTout']) for k in range(1,outputStepping['nDTout']+1)]
