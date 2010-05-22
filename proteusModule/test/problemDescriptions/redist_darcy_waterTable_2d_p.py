from pyadh import *
from pyadh.default_p import *
from math import *

nd = 2

coefficients = RedistanceLevelSet(epsFact=1.5,nModelId=1,rdModelId=-1)
#now define the Dirichlet boundary conditions

def getDBC(x):
    pass
    
dirichletConditions = {0:getDBC}

weakDirichletConditions = {0:RedistanceLevelSet.setZeroLSweakDirichletBCs}
#weakDirichletConditions = None


initialConditions  = None

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 0.75e4
