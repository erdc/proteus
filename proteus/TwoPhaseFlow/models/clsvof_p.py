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
myparams = params.Models.clsvof # model parameters

# MESH #
meshparams = params.mesh
genMesh = meshparams.genMesh

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
                                   useMetrics=myparams.useMetrics,
                                   epsFactHeaviside=myparams.epsFactHeaviside,
                                   epsFactDirac=myparams.epsFactDirac,
                                   epsFactRedist=myparams.epsFactRedist,
                                   lambdaFact=myparams.lambdaFact,
                                   outputQuantDOFs=myparams.outputQuantDOFs,
                                   computeMetrics=myparams.computeMetrics)
coefficients.variableNames=['phi']
name="clsvof"

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions  = {0: initialConditions['clsvof']}

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
