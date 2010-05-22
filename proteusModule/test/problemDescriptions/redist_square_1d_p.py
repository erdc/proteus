from pyadh import *
from pyadh.default_p import *
from math import *
from square import *
from pyadh import RDLS
if useRDLS:
    LevelModelType = RDLS.OneLevelRDLS
class SquareWave1D:
    def __init__(self):
        self.radius = 0.25
        self.velocity=1.0
    def uOfXT(self,x,t):
        center = (0.5 + t)%1.0
        return self.radius - abs(x[0] - center)

analyticalSolution = {0:SquareWave1D()}

coefficients = RedistanceLevelSet(applyRedistancing=applyRedistancing,
                                  epsFact=epsFactRedistance,
                                  nModelId=0,
                                  rdModelId=1)

#now define the Dirichlet boundary conditions

def getDBC(x):
    pass
    
dirichletConditions = {0:getDBC}

def getPDBC(x):
    if x[0] == 0.0 or x[0] == 1.0:
        return numpy.array([0.0,0.0,0.0])

#periodic boundary conditions

periodicDirichletConditions = {0:getPDBC}
if LevelModelType == RDLS.OneLevelRDLS:# or True:
    weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCs}
    periodicDirichletConditions = None
    #weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}
else:
    weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}
    

initialConditions  = {0:analyticalSolution[0]}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
