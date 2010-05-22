from pyadh import *
from pyadh.default_p import *
from math import *
from twp_step3d import *
"""
The redistancing equation in the step test problem.
"""
##

##\ingroup test
#\brief The redistancing equation in the step test problem.
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
    weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}
#weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs2}

class Flat_phi:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        return x[2] - self.waterLevel

initialConditions  = {0:Flat_phi(waterLevel)}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
