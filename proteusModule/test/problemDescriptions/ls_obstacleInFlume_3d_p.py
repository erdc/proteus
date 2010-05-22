from pyadh import *
from pyadh.default_p import *
from obstacleInFlume3d import *
from pyadh import NCLS

LevelModelType = NCLS.OneLevelNCLS

if applyCorrection:
    coefficients = NCLevelSetCoefficients(V_model=0,RD_model=3,ME_model=1)
elif applyRedistancing:
    coefficients = NCLevelSetCoefficients(V_model=0,RD_model=2,ME_model=1)
else:
    coefficients = NCLevelSetCoefficients(V_model=0,ME_model=1)

def getDBC_ls(x,flag):
    pass

dirichletConditions = {0:getDBC_ls}


class Shock_phi:
    def uOfXT(self,x,t):
        return shockSignedDistance(x)

initialConditions  = {0:Shock_phi()}
    
fluxBoundaryConditions = {0:'outFlow'}
    
advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
