from proteus import *
from proteus.default_p import *
from math import *
from risingBubble import *
from proteus.mprans import RDLS3P
"""
The redistancing equation in the sloshbox test problem.
"""

LevelModelType = RDLS3P.LevelModel
coefficients = RDLS3P.Coefficients(applyRedistancing=applyRedistancing,
                                   epsFact=epsFact_redistance,
                                   nModelId=LS_model,
                                   rdModelId=RD_model,
                                   useMetrics=useMetrics,
                                   ELLIPTIC_REDISTANCING=ELLIPTIC_REDISTANCING,
                                   ELLIPTIC_REDISTANCING_TYPE=ELLIPTIC_REDISTANCING_TYPE,
                                   alpha=alpha_REDISTANCING)

def getDBC_rd(x,flag):
    pass

dirichletConditions     = {0:getDBC_rd}
weakDirichletConditions = {0:RDLS3P.setZeroLSweakDirichletBCsSimple}

advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}

class PHI_IC:
    def uOfXT(self,x,t):
        return signedDistanceToBubble(x)

initialConditions  = {0:PHI_IC()}
