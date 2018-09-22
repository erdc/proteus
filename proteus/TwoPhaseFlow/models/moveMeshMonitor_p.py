from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus.mprans import MoveMeshMonitor
import numpy as np
from proteus import Context

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()

myTpFlowProblem = ct.myTpFlowProblem 
initialConditions   = myTpFlowProblem.initialConditions
boundaryConditions  = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd

# DOMAIN #
domain = myTpFlowProblem.domain

params = myTpFlowProblem.Parameters
mparams = myTpFlowProblem.Models # model parameters
pparams = myTpFlowProblem.physical # physical parameters

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
func = mparams.movemeshmonitor['func']
he_max = mparams.movemeshmonitor['he_max']
he_max = mparams.movemeshmonitor['he_min']
ntimes_solved = mparams.movemeshmonitor['ntimes_solved']
fixedNodeMaterialTypes = mparams.movemeshmonitor['fixedNodeMaterialTypes']
fixedElementMaterialTypes = mparams.movemeshmonitor['fixedElementMaterialTypes']
noNodeVelocityMaterialTypes = mparams.movemeshmonitor['noNodeVelocityMaterialTypes']
nSmoothOut = mparams.movemeshmonitor['nSmoothOut']
epsFact = mparams.movemeshmonitor['epsFact']
epsTimeStep = mparams.movemeshmonitor['epsTimeStep']
grading = mparams.movemeshmonitor['grading']
grading_type = mparams.movemeshmonitor['grading_type']
resetNodeVelocityArray = mparams.movemeshmonitor['resetNodeVelocityArray']
useLS = mparams.movemeshmonitor['useLS']
do_firstStep = mparams.movemeshmonitor['do_firstStep']
func = mparams.movemeshmonitor['func']

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
ME_model = mparams.movemeshmonitor['index']
assert ME_model is not None, 'vof model index was not set!'
if mparams.movemeshmonitor['useLS'] is True:
    LS_MODEL = mparams.ncls['index']
else:
    LS_MODEL = None
initialConditions = None

analyticalSolution = {}

coefficients = MoveMeshMonitor.Coefficients(func=func,
                                            nd=nd,
                                            he_max=he_max,
                                            he_min=he_min,
                                            ntimes_solved=movemesh_ntimes,
                                            LS_MODEL=LS_MODEL,
                                            ME_MODEL=ME_MODEL,
                                            fixedNodeMaterialTypes=fixedNodeMaterialTypes,
                                            fixedElementMaterialTypes=fixedElementMaterialTypes,
                                            noNodeVelocityNodeMaterialTypes=noNodeVelocityNodeMaterialTypes,
                                            nSmoothOut=nSmoothOut,
                                            epsTimeStep=epsTimeStep,
                                            epsFact_density=epsFact,
                                            grading=grading,
                                            grading_type=grading_type,
                                            scale_with_nd=scale_with_nd,
                                            do_firstStep=do_firstStep)

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions = None

# ***************************************** #
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
def getDBC(x, flag):
    if flag != 0:
        return None
dirichletConditions = {0: getDBC}

def getDFBC(x, flag):
    return lambda x, t: 0.
diffusiveFluxBoundaryConditions = {0: {0: getDFBC}}
