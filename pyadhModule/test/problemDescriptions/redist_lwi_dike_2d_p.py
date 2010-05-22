from pyadh import *
from pyadh.default_p import *
from math import *
from lwi_dike import *

coefficients = RedistanceLevelSet(applyRedistancing=applyRedistancing,
                                  epsFact=epsFact_redistance,
                                  nModelId=2,
                                  rdModelId=4,
                                  vofModelId=3,
                                  massCorrModelId=5)

#now define the Dirichlet boundary conditions

def getDBC(x):
    pass
    
dirichletConditions = {0:getDBC}

#weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}

initialConditions  = None

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
