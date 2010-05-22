from pyadh import *
from pyadh.default_p import *
from math import *
from embankmentRunup import *
"""
The redistancing equation in the sloshbox test problem.
"""
##

##\ingroup test
#\brief The redistancing equation in the sloshbox test problem.
#
if applyCorrection:
    coefficients = RedistanceLevelSet(applyRedistancing=applyRedistancing,
                                      epsFact=epsFact_redistance,
                                      nModelId=2,
                                      rdModelId=4)
else:
    coefficients = RedistanceLevelSet(applyRedistancing=applyRedistancing,
                                      epsFact=epsFact_redistance,
                                      nModelId=2,
                                      rdModelId=3)

#now define the Dirichlet boundary conditions

def getDBC_rd(x,flag):
    pass
    
dirichletConditions = {0:getDBC_rd}
if rd_freezeLS:
    weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}
    #weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs2}

initialConditions  = None

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
