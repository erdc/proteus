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
movingDomain = myTpFlowProblem.movingDomain

params = myTpFlowProblem.Parameters
mparams = params.Models # model parameters
pparams = params.physical # physical parameters

# MESH #
meshparams = params.mesh
genMesh = meshparams.genMesh

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
myparams = mparams.moveMeshMonitor
func = myparams.func
he_max = myparams.he_max
he_min = myparams.he_min
ntimes_solved = myparams.ntimes_solved
fixedNodeMaterialTypes = myparams.fixedNodeMaterialTypes
fixedElementMaterialTypes = myparams.fixedElementMaterialTypes
noNodeVelocityNodeMaterialTypes = myparams.noNodeVelocityNodeMaterialTypes
nSmoothOut = myparams.nSmoothOut
epsFact = myparams.epsFact
epsTimeStep = myparams.epsTimeStep
grading = myparams.grading
grading_type = myparams.grading_type
resetNodeVelocityArray = myparams.resetNodeVelocityArray
useLS = myparams.useLS
do_firstStep = myparams.do_firstStep
scale_with_nd = myparams.scale_with_nd

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
ME_MODEL = mparams.moveMeshMonitor['index']
assert ME_MODEL is not None, 'moveMeshMonitor model index was not set!'
if mparams.moveMeshMonitor['useLS'] is True:
    LS_MODEL = mparams.ncls['index']
else:
    LS_MODEL = None
initialConditions = None

analyticalSolution = {}

coefficients = MoveMeshMonitor.Coefficients(func=func,
                                            nd=nd,
                                            he_max=he_max,
                                            he_min=he_min,
                                            ntimes_solved=ntimes_solved,
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
