from math import *
from proteus import *
from proteus.default_p import *
from NS_convergence import *

#domain = ctx.domain
#nd = ctx.nd
name = "pressureincrement"

from proteus.mprans import PresInc
coefficients=PresInc.Coefficients(rho_f_min = (1.0-1.0e-8)*rho_1,
                                  rho_s_min = (1.0-1.0e-8)*rho_s,
                                  nd = nd,
                                  modelIndex=PINC_model,
                                  fluidModelIndex=V_model,
                                  fixNullSpace=fixNullSpace_PresInc, 
                                  INTEGRATE_BY_PARTS_DIV_U=INTEGRATE_BY_PARTS_DIV_U_PresInc)
LevelModelType = PresInc.LevelModel

#Always set to None for now
def getDBC_phi(x,flag):    
    None
    
#We should not add any advective flux. The quantity u.n should be determined by the momentum equation
def getAdvectiveFlux_qt(x,flag):
    None

def getDiffusiveFlux_phi(x,flag):
    return lambda x,t: 0.

class getIBC_phi:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:getIBC_phi()}
dirichletConditions = {0:getDBC_phi}
advectiveFluxBoundaryConditions = {0:getAdvectiveFlux_qt}
diffusiveFluxBoundaryConditions = {0:{0:getDiffusiveFlux_phi}}
