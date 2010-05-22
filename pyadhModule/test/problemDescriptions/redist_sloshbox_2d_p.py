from pyadh import *
from pyadh.default_p import *
from math import *
from sloshbox import *
from pyadh import RDLS
"""
The redistancing equation in the sloshbox test problem.
"""
##
#if useRDLS:
#    LevelModelType = RDLS.OneLevelRDLS


##\ingroup test
#\brief The redistancing equation in the sloshbox test problem.
#
if applyCorrection:
    coefficients = RedistanceLevelSet(applyRedistancing=applyRedistancing,
                                      epsFact=epsFact_redistance,
                                      nModelId=1,
                                      rdModelId=3)
else:
    coefficients = RedistanceLevelSet(applyRedistancing=applyRedistancing,
                                      epsFact=epsFact_redistance,
                                      nModelId=1,
                                      rdModelId=2)

#now define the Dirichlet boundary conditions

def getDBC_rd(x,flag):
    pass
    
dirichletConditions = {0:getDBC_rd}

if freezeLevelSet:
    if LevelModelType == RDLS.OneLevelRDLS:
        weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCs}
    else:
        weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}
#weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs2}
#weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs3}

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

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
