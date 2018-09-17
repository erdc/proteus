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
params = ct.params


params.initializeParameters()

# **************************** #
# ********** pnList ********** #
# **************************** #
# List of p/n files
pnList = [None for i in range(ct.params.nModels)]
for i in range(ct.params.nModels):
    model = ct.params.models[i]
    print(model['name'], i, ct.params.nModels)
    pnList[model['index']] = (model['name']+'_p', model['name']+'_n')
    
# ****************************************** #
# ********** TIME STEP CONTROLLER ********** #
# ****************************************** #
systemStepControllerType = Sequential_MinAdaptiveModelStep
if ct.opts.dt_fixed:
    systemStepControllerType = Sequential_FixedStep
    dt_system_fixed = ct.opts.dt_fixed
if params.rans3p['index'] is not None: #rans3p
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
tnList=[0.,ct.outputStepping['dt_init']]+[float(k)*ct.outputStepping['final_time']/float(ct.outputStepping['nDTout']) for k in range(1,ct.outputStepping['nDTout']+1)]
outputStepping = ct.outputStepping
if outputStepping['dt_output'] is None:
    if outputStepping['dt_fixed'] > 0:
        archiveFlag = ArchiveFlags.EVERY_USER_STEP
        if outputStepping['dt_init'] < outputStepping['dt_fixed']:
            tnList = [0., outputStepping['dt_init'], outputStepping['dt_fixed'], outputStepping['final_time']]
        else:
            tnList = [0., outputStepping['dt_fixed'], outputStepping['final_time']]
    else:
          tnList = [0., outputStepping['dt_init'], outputStepping['final_time']]
systemStepExact = False
archiveFlag = ArchiveFlags.EVERY_USER_STEP
if ct.opts.archiveAllSteps is True:
    archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
