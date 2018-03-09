from proteus import *
from proteus.default_p import *
from math import *
from vortex2D import *
from proteus.mprans import RDLS3P
import ncls3P_p
name = soname+"_rdls"
LevelModelType = RDLS3P.LevelModel

coefficients = RDLS3P.Coefficients(applyRedistancing=True,
                                   epsFact=epsFactRedistance,
                                   nModelId=0,
                                   rdModelId=1,
                                   useMetrics=useMetrics,
                                   ELLIPTIC_REDISTANCING=ct.ELLIPTIC_REDISTANCING,
                                   ELLIPTIC_REDISTANCING_TYPE=ct.ELLIPTIC_REDISTANCING_TYPE,
                                   alpha=ct.alpha)

#now define the Dirichlet boundary conditions
def getDBC(x,flag):
    pass

dirichletConditions = {0:getDBC}

if LevelModelType == RDLS3P.LevelModel:
    weakDirichletConditions = {0:RDLS3P.setZeroLSweakDirichletBCsSimple}
else:
    weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}

#weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs2}
#weakDirichletConditions = None

initialConditions  = ncls3P_p.initialConditions
fluxBoundaryConditions = {0:'noFlow'}
advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}
