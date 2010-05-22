from pyadh import *
from pyadh.default_p import *
from threep_cylinder_md_3d import *
from pyadh import NCLS

LevelModelType = NCLS.OneLevelNCLS

coefficients = NCLevelSetCoefficients(V_model=0,RD_model=3,ME_model=1)

def getDBC_ls(x,flag):
    pass

dirichletConditions = {0:getDBC_ls}

class PerturbedSurface_phi:
    def __init__(self,waterLevel,slopeAngle):
        self.waterLevel=waterLevel
        self.slopeAngle=slopeAngle
    def uOfXT(self,x,t):
        signedDistance = x[2] - self.waterLevel
        return signedDistance

initialConditions  = {0:PerturbedSurface_phi(waterLevel,0.0)}
    
fluxBoundaryConditions = {0:'outFlow'}
    
advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
