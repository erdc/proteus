from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus import Context
from proteus.mprans import CLSVOF

# ***************************** #
# ********** CONTEXT ********** #
# ***************************** #
ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions

# ******************************** #
# ********** PARAMETERS ********** #
# ******************************** #
useMetrics = ct.clsvof_parameters['useMetrics']
epsFactHeaviside = ct.clsvof_parameters['epsFactHeaviside']
epsFactDiract = ct.clsvof_parameters['epsFactDirac']
epsFactRedist = ct.clsvof_parameters['epsFactRedist']
lambdaFact = ct.clsvof_parameters['lambdaFact']
outputQuantDOFs = ct.clsvof_parameters['outputQuantDOFs']
computeMetrics = ct.clsvof_parameters['computeMetrics']

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
if ct.opts.ns_model==0: #rans2p
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
                                   computeMetrics=computeMetrics)
coefficients.variableNames=['phi']
name="clsvof"

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions  = {0:ct.clsvof_init_cond()}

# ***************************************** #    
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
dirichletConditions = {0: ct.clsvof_DBC}
advectiveFluxBoundaryConditions = {0: ct.clsvof_AFBC}
diffusiveFluxBoundaryConditions = {0:{0: ct.clsvof_DFBC}}
