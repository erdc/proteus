from __future__ import absolute_import
from proteus import *
from proteus.default_p import *
from math import *
try:
    from .vortex2D import *
    from . import ls_vortex_2d_p
except:
    from vortex2D import *
    import ls_vortex_2d_p
    
from proteus.mprans import RDLS
name = soname+"_rdls"
LevelModelType = RDLS.LevelModel

coefficients = RDLS.Coefficients(applyRedistancing=applyRedistancing,
                                 epsFact=epsFactRedistance,
                                 nModelId=0,
                                 rdModelId=1,
                                 useMetrics=useMetrics,
                                 useExact=useExact)

#now define the Dirichlet boundary conditions

def getDBC(x,flag):
    pass

dirichletConditions = {0:getDBC}

if LevelModelType == RDLS.LevelModel:
    #weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCs}
    weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCsSimple}
else:
    weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}

#weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs2}
#weakDirichletConditions = None


initialConditions  = ls_vortex_2d_p.initialConditions

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
