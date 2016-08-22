from proteus import *
from proteus.default_p import *
from math import *
from floating_bar import *
from proteus.mprans import RDLS
from proteus import Context
ct = Context.get()

"""
The redistancing equation in the sloshbox test problem.
"""

LevelModelType = RDLS.LevelModel

coefficients = RDLS.Coefficients(applyRedistancing=applyRedistancing,
                                 epsFact=epsFact_redistance,
                                 nModelId=int(ct.movingDomain)+2,
                                 rdModelId=int(ct.movingDomain)+3,
                                 useMetrics=useMetrics,
                                 backgroundDiffusionFactor=ct.backgroundDiffusionFactor)

def getDBC_rd(x,flag):
    pass

dirichletConditions     = {0:getDBC_rd}
weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCsSimple}

advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}

class PHI_IC:
    def uOfXT(self,x,t):
        return x[2]-waterLevel

initialConditions  = {0:PHI_IC()}
