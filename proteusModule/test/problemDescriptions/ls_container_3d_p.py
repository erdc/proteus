from pyadh import *
from pyadh.default_p import *
from container import *
from pyadh import NCLS

LevelModelType = NCLS.OneLevelNCLS

#coefficients = NCLevelSetCoefficients(V_model=0,RD_model=3,ME_model=1)
coefficients = NCLevelSetCoefficients(V_model=0,ME_model=1)

def getDBC_ls(x,flag):
    pass

dirichletConditions = {0:getDBC_ls}

class Flat_phi:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        signedDistance = x[2] - self.waterLevel
        return signedDistance

initialConditions  = {0:Flat_phi(waterLevel)}
    
fluxBoundaryConditions = {0:'outFlow'}
    
advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
