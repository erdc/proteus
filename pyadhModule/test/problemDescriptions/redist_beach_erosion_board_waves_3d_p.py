from pyadh import *
from pyadh.default_p import *
from math import *
from beach_erosion_board_waves_3d import *
from pyadh import RDLS
if useRDLS:
    LevelModelType = RDLS.OneLevelRDLS
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
if rd_freezeLS:
    if LevelModelType == RDLS.OneLevelRDLS:
        weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCs}
    else:
        weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}

initialConditions  = None

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
