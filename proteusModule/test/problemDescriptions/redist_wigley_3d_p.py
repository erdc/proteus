from pyadh import *
from pyadh.default_p import *
from math import *
from wigley import *
from pyadh import RDLS

LevelModelType = RDLS.OneLevelRDLS
coefficients = RedistanceLevelSet(applyRedistancing=applyRedistancing,
                                  epsFact=epsFact_redistance,
                                  nModelId=1,
                                  rdModelId=3)

def getDBC_rd(x,flag):
    pass
    
dirichletConditions = {0:getDBC_rd}

if freezeLevelSet:
    if LevelModelType == RDLS.OneLevelRDLS:
        weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCs}
    else:
        weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}
#weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs2}

class Flat_phi:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        signedDistance = x[2] - self.waterLevel
        return signedDistance

initialConditions  = {0:Flat_phi(waterLevel)}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
