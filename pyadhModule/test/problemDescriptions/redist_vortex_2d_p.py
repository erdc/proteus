from pyadh import *
from pyadh.default_p import *
from math import *
from vortex import *
from pyadh import RDLS
if tryRDLS:
    LevelModelType = RDLS.OneLevelRDLS

coefficients = RedistanceLevelSet(applyRedistancing=applyRedistancing,
                                  epsFact=epsFactRedistance,
                                  nModelId=0,
                                  rdModelId=1)

#now define the Dirichlet boundary conditions

def getDBC(x,flag):
    pass
    
dirichletConditions = {0:getDBC}

if LevelModelType == RDLS.OneLevelRDLS:
    #weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCs}
    weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCsSimple}
else:
    #weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCsSimple}
    weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs3}


#weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs2}
#weakDirichletConditions = None


initialConditions  = None

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
