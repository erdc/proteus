"""
Split operator module for multiphase: CLSVOF with RANS2P
"""
from __future__ import absolute_import
from builtins import range
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

# ***************************** #
# ********** CONTEXT ********** #
# ***************************** #
# name = "TwoPhaseFlow"
case = __import__(name)
Context.setFromModule(case)
ct = Context.get()
ct.myTpFlowProblem.initializeAll()

params = ct.myTpFlowProblem.Parameters
params.initializeParameters()

# READ FROM myTpFlowProblem #
assert hasattr(ct,'myTpFlowProblem'), "Create myTpFlowProblem from TwoPhaseFlowProblem"
ns_model = ct.myTpFlowProblem.ns_model
outputStepping = ct.myTpFlowProblem.outputStepping

# **************************** #
# ********** pnList ********** #
# **************************** #
# List of p/n files
pnList = [None for i in range(params.nModels)]
for i in range(params.nModels):
    model = params.models_list[i]
    print(model['name'], i, params.nModels)
    pnList[model['index']] = (model['name']+'_p', model['name']+'_n')

# ****************************************** #
# ********** TIME STEP CONTROLLER ********** #
# ****************************************** #
systemStepControllerType = Sequential_MinAdaptiveModelStep
if ct.myTpFlowProblem.outputStepping['dt_fixed']:
    systemStepControllerType = Sequential_FixedStep
    dt_system_fixed = ct.opts.dt_fixed
if params.Models.rans3p['index'] is not None: #rans3p
    PINIT_model = params.Models.pressureInitial['index']
    assert PINIT_model is not None, 'must set pressureInitial model index when using rans3p'
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
systemStepExact = False
archiveFlag = ArchiveFlags.EVERY_USER_STEP
# if ct.opts.archiveAllSteps is True:
#     archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP