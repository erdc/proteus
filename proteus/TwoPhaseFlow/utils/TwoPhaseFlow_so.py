"""
Split operator module for multiphase: CLSVOF with RANS2P
"""
from __future__ import absolute_import
from builtins import range
import os
from proteus.default_so import *
from proteus import Context
import sys

<<<<<<< HEAD
# ***************************** #
# ********** CONTEXT ********** #
# ***************************** #
name = "TwoPhaseFlow"
pathMyTpFlowProblem = None
=======
pathMyTpFlowProblem = None
name = 'TwoPhaseFlow'

>>>>>>> TwoPhaseFlow

for i in range(len(sys.argv)):
    if '-f' == sys.argv[i]:
        assert sys.argv[i+1][-3:], "fileName must end with .py"
        name = sys.argv[i+1][:-3]
    if '--path' == sys.argv[i]:
        pathMyTpFlowProblem = sys.argv[i+1]

# Load module
if pathMyTpFlowProblem is not None:
    sys.path.append(pathMyTpFlowProblem)
<<<<<<< HEAD
=======

# ***************************** #
# ********** CONTEXT ********** #
# ***************************** #
# name = "TwoPhaseFlow"
>>>>>>> TwoPhaseFlow
case = __import__(name)

# Create context
Context.setFromModule(case)
ct = Context.get()
<<<<<<< HEAD
=======
ct.myTpFlowProblem.initializeAll()

params = ct.myTpFlowProblem.Parameters
>>>>>>> TwoPhaseFlow


# READ FROM myTpFlowProblem #
assert hasattr(ct,'myTpFlowProblem'), "Create myTpFlowProblem from TwoPhaseFlowProblem"
ns_model = ct.myTpFlowProblem.ns_model
outputStepping = ct.myTpFlowProblem.outputStepping

# **************************** #
# ********** pnList ********** #
# **************************** #
<<<<<<< HEAD
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
=======
# List of p/n files
pnList = [None for i in range(params.nModels)]
for i in range(params.nModels):
    model = params.models_list[i]
    print(model['name'], i, params.nModels)
    pnList[model['index']] = (model['name']+'_p', model['name']+'_n')
>>>>>>> TwoPhaseFlow

# ****************************************** #
# ********** TIME STEP CONTROLLER ********** #
# ****************************************** #
systemStepControllerType = Sequential_MinAdaptiveModelStep
<<<<<<< HEAD
if ns_model==1: #rans3p
    PINIT_model=4
=======
if ct.myTpFlowProblem.outputStepping['dt_fixed']:
    systemStepControllerType = Sequential_FixedStep
    dt_system_fixed = ct.opts.dt_fixed
if params.Models.rans3p['index'] is not None: #rans3p
    PINIT_model = params.Models.pressureInitial['index']
    assert PINIT_model is not None, 'must set pressureInitial model index when using rans3p'
>>>>>>> TwoPhaseFlow
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
<<<<<<< HEAD
tnList=[0.,outputStepping['dt_init']]+[float(k)*outputStepping['final_time']/float(outputStepping['nDTout']) for k in range(1,outputStepping['nDTout']+1)]
=======
outputStepping = ct.myTpFlowProblem.outputStepping
tnList=[0.,outputStepping['dt_init']]+[float(k)*outputStepping['final_time']/float(outputStepping['nDTout']) for k in range(1,outputStepping['nDTout']+1)]
if outputStepping['dt_output'] is None:
    if outputStepping['dt_fixed'] > 0:
        if outputStepping['dt_init'] < outputStepping['dt_fixed']:
            tnList = [0., outputStepping['dt_init'], outputStepping['dt_fixed'], outputStepping['final_time']]
        else:
            tnList = [0., outputStepping['dt_fixed'], outputStepping['final_time']]
    else:
        tnList = [0., outputStepping['dt_init'], outputStepping['final_time']]
systemStepExact = ct.myTpFlowProblem.outputStepping.systemStepExact

archiveFlag = ArchiveFlags.EVERY_USER_STEP
if ct.myTpFlowProblem.archiveAllSteps is True:
    archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
>>>>>>> TwoPhaseFlow
