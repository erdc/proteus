from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from math import *
from dambreak import *
from proteus.mprans import CLSVOF

LevelModelType = CLSVOF.LevelModel
#coefficients = MyCoefficients(
coefficients = CLSVOF.Coefficients(
    V_model=V_model,
    useMetrics=useMetrics,
    timeOrder=timeOrder,
    epsFactHeaviside=1.5,
    epsFactDirac=1.5,
    lambdaFact=1.0,
    outputQuantDOFs=True,
    computeMetrics=0)
coefficients.variableNames=['phi']

#####################
# INITIAL CONDITION #
#####################
class init_cond:
    def uOfXT(self,x,t):
        return signedDistanceToBubble(x)
initialConditions  = {0:init_cond()}

#######################
# BOUNDARY CONDITIONS #
#######################
def getDBC_vof(x,flag):
    if flag == boundaryTags['top'] and openTop:
        return lambda x,t: 1.0
dirichletConditions = {0:getDBC_vof}

def getAFBC_vof(x,flag):
    if flag != boundaryTags['top'] or not openTop:
        return lambda x,t: 0.0
advectiveFluxBoundaryConditions = {0:getAFBC_vof}
diffusiveFluxBoundaryConditions = {0:{}}
