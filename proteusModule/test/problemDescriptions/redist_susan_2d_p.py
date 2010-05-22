from pyadh import *
from pyadh.default_p import *
from math import *
from susan import *
"""
The redistancing equation used by the SUSAN MAERSK test problem
"""
##

##\ingroup test
#\brief The redistancing equation used by the SUSAN MAERSK test problem
coefficients = RedistanceLevelSet(epsFact=0.5*epsFact,nModelId=1,rdModelId=2)

#now define the Dirichlet boundary conditions

def getDBC(x):
    if x[0] >= (2.0e2-1.0e-8):
        return lambda x,t: x[1] - waterLevel

dirichletConditions = {0:getDBC}

weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}
#weakDirichletConditions = None

initialConditions  = None

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
