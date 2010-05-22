from pyadh import *
from pyadh.default_p import *
from obstacleInFlume3d import *
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

def getDBC_p_obstacleInFlume(x,flag):
    if flag == boundaryTags['right']:
        return lambda x,t: -coefficients.g[2]*(rho_0*(height - x[2])
                                               -(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,height-waterLevel)
                                               +(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,x[2]-waterLevel))
def getDBC_u_obstacleInFlume(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: Um
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0

def getDBC_v_obstacleInFlume(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0
def getDBC_w_obstacleInFlume(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p_obstacleInFlume,
                       1:getDBC_u_obstacleInFlume,
                       2:getDBC_v_obstacleInFlume,
                       3:getDBC_w_obstacleInFlume}

def getAFBC_p_obstacleInFlume(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: -Um#*velRamp(t)
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0
    if flag == boundaryTags['front']:
        return lambda x,t: 0.0
    if flag == boundaryTags['back']:
        return lambda x,t: 0.0

def getAFBC_u_obstacleInFlume(x,flag):
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0
    if flag == boundaryTags['front']:
        return lambda x,t: 0.0
    if flag == boundaryTags['back']:
        return lambda x,t: 0.0
#    if flag == boundaryTags['obstacle']:
#        return lambda x,t: 0.0

def getAFBC_v_obstacleInFlume(x,flag):
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0
    if flag == boundaryTags['front']:
        return lambda x,t: 0.0
    if flag == boundaryTags['back']:
        return lambda x,t: 0.0
#    if flag == boundaryTags['obstacle']:
#        return lambda x,t: 0.0

def getAFBC_w_obstacleInFlume(x,flag):
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0
    if flag == boundaryTags['front']:
        return lambda x,t: 0.0
    if flag == boundaryTags['back']:
        return lambda x,t: 0.0
#    if flag == boundaryTags['obstacle']:
#        return lambda x,t: 0.0

def getDFBC_u_obstacleInFlume(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    if flag == boundaryTags['right']:
        return lambda x,t: 0.0
    if flag == boundaryTags['front']:
        return lambda x,t: 0.0
    if flag == boundaryTags['back']:
        return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0
#    if flag == boundaryTags['obstacle']:
#        return lambda x,t: 0.0

def getDFBC_v_obstacleInFlume(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    if flag == boundaryTags['right']:
        return lambda x,t: 0.0
    if flag == boundaryTags['front']:
        return lambda x,t: 0.0
    if flag == boundaryTags['back']:
        return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0
 #   if flag == boundaryTags['obstacle']:
 #       return lambda x,t: 0.0

def getDFBC_w_obstacleInFlume(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0
    if flag == boundaryTags['right']:
        return lambda x,t: 0.0
    if flag == boundaryTags['front']:
        return lambda x,t: 0.0
    if flag == boundaryTags['back']:
        return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0
  #  if flag == boundaryTags['obstacle']:
  #      return lambda x,t: 0.0

fluxBoundaryConditions = {0:'mixedFlow',
                          1:'mixedFlow',
                          2:'mixedFlow',
                          3:'mixedFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_obstacleInFlume,
                                    1:getAFBC_u_obstacleInFlume,
                                    2:getAFBC_v_obstacleInFlume,
                                    3:getAFBC_w_obstacleInFlume}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u_obstacleInFlume},
                                   2:{2:getDFBC_v_obstacleInFlume},
                                   3:{3:getDFBC_w_obstacleInFlume}}

class Steady_p:
    def uOfXT(self,x,t):
        return -coefficients.g[2]*(rho_0*(height - x[2])
                                   -(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,height-waterLevel)
                                   +(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,x[2]-waterLevel))

class Steady_u:
    def uOfXT(self,x,t):
        return Um

class Steady_v:
    def uOfXT(self,x,t):
        return 0.0

class Steady_w:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v(),
                     3:Steady_w()}
