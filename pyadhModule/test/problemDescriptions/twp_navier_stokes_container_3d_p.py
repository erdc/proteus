"""
Incompressible Navier-Stokes flow around a square obstacle in 2D.
"""
from pyadh import *
from pyadh.default_p import *
import sys
from math import *
from container import *
from pyadh.ctransportCoefficients import smoothedHeaviside
from pyadh.ctransportCoefficients import smoothedHeaviside_integral
from pyadh import RANS2P

LevelModelType = RANS2P.OneLevelRANS2P
coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=epsFact_viscosity,
                                             sigma=sigma_01,
                                             rho_0=rho_0,
                                             nu_0=nu_0,
                                             rho_1=rho_1,
                                             nu_1=nu_1,
                                             g=g,
                                             nd=nd,
                                             LS_model=1,
                                             epsFact_density=epsFact_density,
                                             stokes=useStokes,
                                             movingDomain=movingDomain)
coefficients.waterLevel=waterLevel

def velRamp(t):
    if t < 0.1*residence_time:
        return 1.0-exp(-25.0*t/(0.1*residence_time))
    else:
        return 1.0

def getDBC_p(x,flag):
    if flag == boundaryTags['right']:
        return lambda x,t: -coefficients.g[2]*(rho_0*(height - x[2])
                                               -(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,height-waterLevel)
                                               +(rho_0-rho_1)*smoothedHeaviside_integral(epsFact_density*he,x[2]-waterLevel))

walls = [boundaryTags['front'],boundaryTags['back'],boundaryTags['bottom'],boundaryTags['top']]
def getDBC_u(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: Um
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0

def getDBC_v(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0

def getDBC_w(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v,
                       3:getDBC_w}

def getAFBC_p(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: -Um
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0
    if flag in walls:
        return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0
def getAFBC_u(x,flag):
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0
    if flag in walls:
        return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0

def getAFBC_v(x,flag):
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0
    if flag in walls:
        return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0
    
def getAFBC_w(x,flag):
    if flag == boundaryTags['obstacle']:
        return lambda x,t: 0.0
    if flag in walls:
        return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0
    
def getDFBC_u(x,flag):
    if flag in walls:
        return lambda x,t: 0.0
    if flag == boundaryTags['right']:
        return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0

def getDFBC_v(x,flag):
    if flag in walls:
        return lambda x,t: 0.0
    if flag == boundaryTags['right']:
        return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0

def getDFBC_w(x,flag):
    if flag in walls:
        return lambda x,t: 0.0
    if flag == boundaryTags['right']:
        return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0

fluxBoundaryConditions = {0:'mixedFlow',
                          1:'mixedFlow',
                          2:'mixedFlow',
                          3:'mixedFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getAFBC_u,
                                    2:getAFBC_v,
                                    3:getAFBC_w}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v},
                                   3:{3:getDFBC_w}}

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
