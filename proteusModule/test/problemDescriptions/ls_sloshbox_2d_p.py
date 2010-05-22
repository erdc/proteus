from pyadh import *
from pyadh.default_p import *
from sloshbox import *
"""
The non-conservative level set description of the free surface of a sloshing two-phase flow in a closed box.
"""
## 

##\ingroup test
#\file ls_so_sloshbox_2d_p.py
#
#\brief The non-conservative level set description of the free surface of a sloshing two-phase flow in a closed box.
#
#\todo finish ls_so_sloshbox_2d_p.py doc
analyticalSolutions = None

if applyCorrection:
    coefficients = NCLevelSetCoefficients(V_model=0,RD_model=3,ME_model=1)
elif applyRedistancing:
    coefficients = NCLevelSetCoefficients(V_model=0,RD_model=2,ME_model=1)
else:
    coefficients = NCLevelSetCoefficients(V_model=0,RD_model=-1,ME_model=1)

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
