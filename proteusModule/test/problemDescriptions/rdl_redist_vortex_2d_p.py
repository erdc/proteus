from pyadh import *
from pyadh.default_p import *
from math import *
from vortex import *

coefficients = RedistanceLevelSet(applyRedistancing=applyRedistancing,
                                  epsFact=0.0,
                                  nModelId=0,
                                  rdModelId=2,
                                  vofModelId=1,
                                  massCorrModelId=3)

#now define the Dirichlet boundary conditions

def getDBC(x):
    pass
    
dirichletConditions = {0:getDBC}

weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}
#weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs2}
#weakDirichletConditions = None


initialConditions  = None

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
