from pyadh import *
from pyadh.default_p import *
from damBreak import *

if applyCorrection:
    coefficients = NCLevelSetCoefficients(V_model=1,RD_model=4,ME_model=2)
elif applyRedistancing:
    coefficients = NCLevelSetCoefficients(V_model=1,RD_model=3,ME_model=2)
else:
    coefficients = NCLevelSetCoefficients(V_model=1,RD_model=-1,ME_model=2)

def getDBC_ls(x,flag):
    pass

dirichletConditions = {0:getDBC_ls}

class PerturbedSurface_phi:
    def __init__(self,waterLevel,slopeAngle):
        self.waterLevel=waterLevel
        self.slopeAngle=slopeAngle
    def uOfXT(self,x,t):
        if useShock:
            return shockSignedDistance(x)
        else:
            surfaceNormal = [-sin(self.slopeAngle),cos(self.slopeAngle)]
            signedDistance = (x[0] - 0.5)*surfaceNormal[0]+(x[1] - self.waterLevel)*surfaceNormal[1]
            return signedDistance

initialConditions  = {0:PerturbedSurface_phi(waterLevel,slopeAngle)}
    
fluxBoundaryConditions = {0:'outFlow'}
    
advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
