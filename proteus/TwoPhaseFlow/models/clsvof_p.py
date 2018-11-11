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
clsvof_parameters   = myTpFlowProblem.clsvof_parameters
initialConditions   = myTpFlowProblem.initialConditions
boundaryConditions  = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd
ns_model = myTpFlowProblem.ns_model

# DOMAIN #
domain = myTpFlowProblem.domain

# ******************************** #
# ********** PARAMETERS ********** #
# ******************************** #
useMetrics = clsvof_parameters['useMetrics']
epsFactHeaviside = clsvof_parameters['epsFactHeaviside']
epsFactDiract = clsvof_parameters['epsFactDirac']
epsFactRedist = clsvof_parameters['epsFactRedist']
lambdaFact = clsvof_parameters['lambdaFact']
outputQuantDOFs = clsvof_parameters['outputQuantDOFs']
computeMetrics = clsvof_parameters['computeMetrics']
disc_ICs = clsvof_parameters['disc_ICs']

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
if ns_model==0: #rans2p
    CLSVOF_model=1
    V_model=0
else:
    CLSVOF_model=0
    V_model=1

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
                                   computeMetrics=computeMetrics,
                                   computeMetricsForBubble=False,
                                   disc_ICs=disc_ICs)
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
