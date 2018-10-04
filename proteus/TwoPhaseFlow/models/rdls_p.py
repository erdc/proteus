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
myparams = params.Models.rdls
pparams = params.physical # physical parameters

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
nModelId = mparams.ncls.index
assert nModelId is not None, 'redist model index was not set!'
rdModelId = mparams.rdls.index

LevelModelType = RDLS.LevelModel
coefficients = RDLS.Coefficients(applyRedistancing=myparams.applyRedistancing,
                                 epsFact=myparams.epsFact,
                                 nModelId=nModelId,
                                 rdModelId=rdModelId,
                                 useMetrics=myparams.useMetrics,
                                 backgroundDiffusionFactor=myparams.backgroundDiffusionFactor)

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions  = {0: initialConditions['rdls']}

# ***************************************** #
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
dirichletConditions     = {0: lambda x, flag: None}
weakDirichletConditions = {0: RDLS.setZeroLSweakDirichletBCsSimple}
advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}


