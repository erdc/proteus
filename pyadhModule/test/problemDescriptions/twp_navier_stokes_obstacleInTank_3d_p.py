from pyadh import *
from pyadh.default_p import *
from obstacleInTank3d import *
from pyadh import RANS2P
from pyadh.ctransportCoefficients import smoothedHeaviside_integral

LevelModelType = RANS2P.OneLevelRANS2P
coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=epsFact_viscosity,
                                             sigma=0.0,
                                             rho_0 = rho_0,
                                             nu_0 = nu_0,
                                             rho_1 = rho_1,
                                             nu_1 = nu_1,
                                             g=g,
                                             nd=nd,
                                             LS_model=1,
                                             epsFact_density=epsFact_density,
                                             stokes=useStokes)

def getDBC_p_obstacleInTank(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0
def getDBC_u_obstacleInTank(x,flag):
    return None
def getDBC_v_obstacleInTank(x,flag):
    return None
def getDBC_w_obstacleInTank(x,flag):
    return None

dirichletConditions = {0:getDBC_p_obstacleInTank,
                       1:getDBC_u_obstacleInTank,
                       2:getDBC_v_obstacleInTank,
                       3:getDBC_w_obstacleInTank}

def getAFBC_p_obstacleInTank(x,flag):
    if(flag == boundaryTags['obstacle'] or
       flag == boundaryTags['left'] or
       flag == boundaryTags['right'] or
       flag == boundaryTags['front'] or
       flag == boundaryTags['back'] or
       flag == boundaryTags['bottom']):
        return lambda x,t: 0.0
def getAFBC_u_obstacleInTank(x,flag):
    if(flag == boundaryTags['obstacle'] or
       flag == boundaryTags['left'] or
       flag == boundaryTags['right'] or
       flag == boundaryTags['front'] or
       flag == boundaryTags['back'] or
       flag == boundaryTags['bottom']):
        return lambda x,t: 0.0
def getAFBC_v_obstacleInTank(x,flag):
    if(flag == boundaryTags['obstacle'] or
       flag == boundaryTags['left'] or
       flag == boundaryTags['right'] or
       flag == boundaryTags['front'] or
       flag == boundaryTags['back'] or
       flag == boundaryTags['bottom']):
        return lambda x,t: 0.0
def getAFBC_w_obstacleInTank(x,flag):
    if(flag == boundaryTags['obstacle'] or
       flag == boundaryTags['left'] or
       flag == boundaryTags['right'] or
       flag == boundaryTags['front'] or
       flag == boundaryTags['back'] or
       flag == boundaryTags['bottom']):
        return lambda x,t: 0.0
def getDFBC_u_obstacleInTank(x,flag):
    return lambda x,t: 0.0
def getDFBC_v_obstacleInTank(x,flag):
    return lambda x,t: 0.0
def getDFBC_w_obstacleInTank(x,flag):
    return lambda x,t: 0.0

fluxBoundaryConditions = {0:'mixedFlow',
                          1:'mixedFlow',
                          2:'mixedFlow',
                          3:'mixedFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_obstacleInTank,
                                    1:getAFBC_u_obstacleInTank,
                                    2:getAFBC_v_obstacleInTank,
                                    3:getAFBC_w_obstacleInTank}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u_obstacleInTank},
                                   2:{2:getDFBC_v_obstacleInTank},
                                   3:{3:getDFBC_w_obstacleInTank}}

class Shock_p:
    def uOfXT(self,x,t):
        if shockSignedDistance(x) < 0:
            return -g[2]*(rho_0*(height - x[2])
                          -(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,height-waterColumn_z)
                          +(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,x[2]-waterColumn_z))
        else:
            return -rho_1*g[2]

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Shock_p(),
                     1:AtRest(),
                     2:AtRest(),
                     3:AtRest()}
