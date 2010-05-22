from pyadh import *
from pyadh.default_p import *
from math import *
from flume import *

coefficients = RedistanceLevelSet(applyRedistancing=applyRedistancing,
                                  epsFact=epsFact_redistance,
                                  nModelId=2,
                                  rdModelId=4)

#now define the Dirichlet boundary conditions

def getDBC(x,flag):
    pass
    
dirichletConditions = {0:getDBC}

weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}

initialConditions  = {0:Flat_phi()}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
