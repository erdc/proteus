from pyadh import *
from pyadh.default_p import *
from math import *
from obstacleInFlume3d import *
from pyadh import RDLS
"""
The redistancing equation in the obstacleInFlume test problem.
"""
##
LevelModelType = RDLS.OneLevelRDLS

##\ingroup test
#\brief The redistancing equation in the obstacleInFlume test problem.
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

class Shock_phi:
    def uOfXT(self,x,t):
        return shockSignedDistance(x)

initialConditions  = {0:Shock_phi()}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
