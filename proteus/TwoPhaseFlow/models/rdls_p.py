from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus.mprans import RDLS
from proteus import Context

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()

myTpFlowProblem = ct.myTpFlowProblem 
initialConditions = myTpFlowProblem.initialConditions
boundaryConditions = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd
movingDomain = myTpFlowProblem.movingDomain


# DOMAIN #
domain = myTpFlowProblem.domain

params = myTpFlowProblem.Parameters
mparams = params.Models # model parameters
pparams = params.physical # physical parameters

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
useMetrics = mparams.rdls['useMetrics']
applyRedistancing = mparams.rdls['applyRedistancing']
epsFact = mparams.rdls['epsFact']
backgroundDiffusionFactor = mparams.rdls['backgroundDiffusionFactor']

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
nModelId = mparams.ncls['index']
assert nModelId != None, 'redist model index was not set!'
rdModelId = mparams.rdls['index']

LevelModelType = RDLS.LevelModel
coefficients = RDLS.Coefficients(applyRedistancing=applyRedistancing,
                                 epsFact=epsFact,
                                 nModelId=nModelId,
                                 rdModelId=rdModelId,
                                 useMetrics=useMetrics,
                                 backgroundDiffusionFactor=backgroundDiffusionFactor)

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions  = {0: initialConditions['rdls']}

# ***************************************** #
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
dirichletConditions     = {0: None}
weakDirichletConditions = {0: RDLS.setZeroLSweakDirichletBCsSimple}
advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}


