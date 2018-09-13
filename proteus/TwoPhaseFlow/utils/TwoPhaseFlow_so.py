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
name = "TwoPhaseFlow"
case = __import__(name)
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
