from pyadh import *
from pyadh.default_p import *
from pyadh.ctransportCoefficients import smoothedHeaviside
from sloshbox3d import *
from pyadh import VOF

LevelModelType = VOF.OneLevelVOF
coefficients = VOFCoefficients(LS_model=1,V_model=0,RD_model=3,ME_model=2,epsFact=epsFact_vof)

class PerturbedSurface_H:
    def __init__(self,waterLevel,slopeAngle):
        self.waterLevel=waterLevel
        self.slopeAngle=slopeAngle
    def uOfXT(self,x,t):
        if useShock:
            return smoothedHeaviside(epsFact_consrv_heaviside*he,shockSignedDistance(x))
        else:
            surfaceNormal = [-sin(self.slopeAngle),cos(self.slopeAngle)]
            signedDistance = (x[0] - 0.5)*surfaceNormal[0]+(x[2] - self.waterLevel)*surfaceNormal[1]
            return smoothedHeaviside(epsFact_consrv_heaviside*he,signedDistance)

analyticalSolutions = None
#closedTop=True
def getDBC_vof(x,flag):
    if not closedTop:
        if regularGrid:
            if x[2] > L[2]-1.0e-8:
                return lambda x,t: 1.0
        else:
            if flag == bt['top']:
                return lambda x,t: 1.0

dirichletConditions = {0:getDBC_vof}

initialConditions  = {0:PerturbedSurface_H(waterLevel,slopeAngle)}

#fluxBoundaryConditions = {0:'noFlow'}
fluxBoundaryConditions = {0:'mixedFlow'}

def getAFBC_vof(x,flag):
    if closedTop:
        return lambda x,t: 0.0
    else:
        if regularGrid:
            if x[2] > L[2]-1.0e-8:
                return None
            else:
                return lambda x,t: 0.0
        else:
            if flag == bt['top']:
                return None
            else:
                return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_vof}

diffusiveFluxBoundaryConditions = {0:{}}
