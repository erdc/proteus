from pyadh import *
from pyadh.default_p import *
from twp_step3d import *

if applyCorrection:
    coefficients = NCLevelSetCoefficients(V_model=0,RD_model=3,ME_model=1)
elif applyRedistancing:
    coefficients = NCLevelSetCoefficients(V_model=0,RD_model=2,ME_model=1)
else:
    coefficients = NCLevelSetCoefficients(V_model=0,ME_model=1)

def getDBC_ls(x,flag):
    pass

dirichletConditions = {0:getDBC_ls}

class Flat_phi:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        return x[2] - self.waterLevel
    
initialConditions  = {0:Flat_phi(waterLevel)}
    
fluxBoundaryConditions = {0:'outFlow'}
    
advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
