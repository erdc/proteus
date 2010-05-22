from pyadh import *
from pyadh.default_p import *
from simpleweir import *
"""
The redistancing equation for the simple weir test problem.
"""
##

##\ingroup test
#\brief The redistancing equation for the simple weir test problem.
#
#The model equations are defined by the RedistanceLevelSet class. The initial/boundary conditions are
#\f{eqnarray*}
#\phi(x,y,0) &=& y - h_w \\
#\phi(0,y,t) &=& y - h_w
#\f}
#where $h_w$ is the height of the free surface
#
analyticalSolutions = None
                                                    
coefficients = RedistanceLevelSet(epsFact=0.5*epsFact,nModelId=1,rdModelId=2)

#now define the Dirichlet boundary conditions

def getDBC(x):
    pass
#     if x[0] == 0.0 or x[0] >= flumeEnd - 1.0e-8:
#         return lambda x,t: x[1]-waterLevel
#     if x[0] == 0.0:
#         return lambda x,t:  x[1]-waterLevel
#    pass

dirichletConditions = {0:getDBC}

weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}
#weakDirichletConditions = None

initialConditions  = None

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
