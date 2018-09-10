from __future__ import absolute_import
from builtins import object
from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from math import *
from proteus.mprans import CLSVOF
from proteus import Context

ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions

print ct
CLSVOF_model=1
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

LevelModelType = CLSVOF.LevelModel
coefficients = CLSVOF.Coefficients(V_model=0,
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

#####################
# INITIAL CONDITION #
#####################
initialConditions  = {0:ct.clsvof_init_cond()}

#######################
# BOUNDARY CONDITIONS #
#######################
dirichletConditions = {0: ct.clsvof_DBC}
advectiveFluxBoundaryConditions = {0: ct.clsvof_AFBC}
diffusiveFluxBoundaryConditions = {0:{0: ct.clsvof_DFBC}}
