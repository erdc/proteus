from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus import Context
from proteus.mprans import CLSVOF

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()
<<<<<<< HEAD
myTpFlowProblem = ct.myTpFlowProblem 
clsvof_parameters   = myTpFlowProblem.clsvof_parameters
initialConditions   = myTpFlowProblem.initialConditions
boundaryConditions  = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd
ns_model = myTpFlowProblem.ns_model
=======

myTpFlowProblem = ct.myTpFlowProblem
initialConditions   = myTpFlowProblem.initialConditions
boundaryConditions  = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd
>>>>>>> TwoPhaseFlow

# DOMAIN #
domain = myTpFlowProblem.domain

<<<<<<< HEAD
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
computeMetricsForBubble = clsvof_parameters['computeMetricsForBubble']
disc_ICs = clsvof_parameters['disc_ICs']
=======
ns_model = myTpFlowProblem.ns_model

params = myTpFlowProblem.Parameters
mparams = params.Models # model parameters
pparams = params.physical # physical parameters
myparams = params.Models.clsvof # model parameters

# MESH #
meshparams = params.mesh
genMesh = meshparams.genMesh
>>>>>>> TwoPhaseFlow

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
<<<<<<< HEAD
if ns_model==0: #rans2p
    CLSVOF_model=1
    V_model=0
else:
    CLSVOF_model=0
    V_model=1
=======
CLSVOF_model = mparams.clsvof['index']
if ns_model==0: #rans2p
    V_model = mparams.rans2p['index']
else:
    V_model = mparams.rans3p['index']
>>>>>>> TwoPhaseFlow

# ********************************** #
# ********** COEFFICIENTS ********** #
# ********************************** #
LevelModelType = CLSVOF.LevelModel
coefficients = CLSVOF.Coefficients(V_model=V_model,
                                   ME_model=CLSVOF_model,
<<<<<<< HEAD
                                   useMetrics=useMetrics,
                                   epsFactHeaviside=epsFactHeaviside,
                                   epsFactDirac=epsFactHeaviside,
                                   epsFactRedist=epsFactRedist,
                                   lambdaFact=lambdaFact,
                                   outputQuantDOFs=outputQuantDOFs,
                                   computeMetrics=computeMetrics,
                                   computeMetricsForBubble=computeMetricsForBubble,
                                   disc_ICs=disc_ICs)
=======
                                   useMetrics=myparams.useMetrics,
                                   epsFactHeaviside=myparams.epsFactHeaviside,
                                   epsFactDirac=myparams.epsFactDirac,
                                   epsFactRedist=myparams.epsFactRedist,
                                   lambdaFact=myparams.lambdaFact,
                                   outputQuantDOFs=myparams.outputQuantDOFs,
                                   computeMetrics=myparams.computeMetrics)
>>>>>>> TwoPhaseFlow
coefficients.variableNames=['phi']
name="clsvof"

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions  = {0: initialConditions['clsvof']}

<<<<<<< HEAD
# ***************************************** #    
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
dirichletConditions = {0: boundaryConditions['clsvof_DBC']}
advectiveFluxBoundaryConditions = {0: boundaryConditions['clsvof_AFBC']}
diffusiveFluxBoundaryConditions = {0:{0: boundaryConditions['clsvof_DFBC']}}
=======
# ***************************************** #
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
if domain.useSpatialTools is False or myTpFlowProblem.useBoundaryConditionsModule is False:
    dirichletConditions = {0: boundaryConditions['clsvof_DBC']}
    advectiveFluxBoundaryConditions = {0: boundaryConditions['clsvof_AFBC']}
    diffusiveFluxBoundaryConditions = {0:{0: boundaryConditions['clsvof_DFBC']}}
else:
    dirichletConditions = {0: lambda x, flag: domain.bc[flag].clsvof_dirichlet.init_cython()}
    advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].clsvof_advective.init_cython()}
    diffusiveFluxBoundaryConditions = {0: {0: lambda x, flag: domain.bc[flag].clsvof_diffusive.init_cython()}}
>>>>>>> TwoPhaseFlow
