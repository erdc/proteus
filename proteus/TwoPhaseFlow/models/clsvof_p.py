from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus import Context
from proteus.mprans import CLSVOF

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

ns_model = myTpFlowProblem.ns_model

params = myTpFlowProblem.Parameters
mparams = params.Models # model parameters
pparams = params.physical # physical parameters

# ******************************** #
# ********** PARAMETERS ********** #
# ******************************** #
useMetrics = mparams.clsvof['useMetrics']
epsFactHeaviside = mparams.clsvof['epsFactHeaviside']
epsFactDiract = mparams.clsvof['epsFactDirac']
epsFactRedist = mparams.clsvof['epsFactRedist']
lambdaFact = mparams.clsvof['lambdaFact']
outputQuantDOFs = mparams.clsvof['outputQuantDOFs']
computeMetrics = mparams.clsvof['computeMetrics']

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
CLSVOF_model = mparams.clsvof['index']
if ns_model==0: #rans2p
    V_model = mparams.rans2p['index']
else:
    V_model = mparams.rans3p['index']

# ********************************** #
# ********** COEFFICIENTS ********** #
# ********************************** #
LevelModelType = CLSVOF.LevelModel
coefficients = CLSVOF.Coefficients(V_model=V_model,
                                   ME_model=CLSVOF_model,
                                   useMetrics=useMetrics,
                                   epsFactHeaviside=epsFactHeaviside,
                                   epsFactDirac=epsFactHeaviside,
                                   epsFactRedist=epsFactRedist,
                                   lambdaFact=lambdaFact,
                                   outputQuantDOFs=outputQuantDOFs,
                                   computeMetrics=computeMetrics)
coefficients.variableNames=['phi']
name="clsvof"

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions  = {0: initialConditions['clsvof']}

# ***************************************** #    
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
dirichletConditions = {0: boundaryConditions['clsvof_DBC']}
advectiveFluxBoundaryConditions = {0: boundaryConditions['clsvof_AFBC']}
diffusiveFluxBoundaryConditions = {0:{0: boundaryConditions['clsvof_DFBC']}}
