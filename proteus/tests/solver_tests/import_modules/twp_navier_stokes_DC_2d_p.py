"""
Incompressible Navier-Stokes flow around a cylinder in 2D.
"""
from proteus import *
from proteus.default_p import *
from proteus import Domain
import sys
from twpDC import *
from proteus.mprans import RANS2P

name = 'mprans_test'

bcsTimeDependent = True
LevelModelType = RANS2P.LevelModel
coefficients = RANS2P.Coefficients(epsFact=epsFact_viscosity,
                                   rho_0=rho,
                                   nu_0=nu,
                                   rho_1=rho,
                                   nu_1=nu,
                                   g=g,
                                   nd=nd,
                                   LS_model=None,
                                   epsFact_density=epsFact_density,
                                   stokes=False,#useStokes,
                                   useRBLES=0.0,
                                   useMetrics=0.0,
                                   forceStrongDirichlet=True)



from math import cos,pi
def velRamp(t):
    return (0.41**-2)*sin(pi*t/8.0)

def top(x,flag):
    if x[1]== x0[1]+L[1]:
        return True
    else:
        return False

def sides(x,flag):
    if (x[0]==x0[0] or x[0]==(x0[0]+L[0])):
        return True
    else:
        return False

def bottom(x,flag):
    if (x[1] == x0[1]):
        return True
    else:
        return False

def getDBC_p(x,flag):
    pass

def getDBC_u(x,flag):
    if sides(x,flag) or bottom(x,flag):
        return lambda x,t: 0.0
    elif top(x,flag):
        return lambda x,t: 1.0

def getDBC_v(x,flag):
    if (top(x,flag) or sides(x,flag) or bottom(x,flag)):
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v}

def getAFBC_p(x,flag):
    pass

def getAFBC_u(x,flag):
    pass

def getAFBC_v(x,flag):
    pass

def getDFBC_u(x,flag):
    pass

def getDFBC_v(x,flag):
    pass

advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getAFBC_u,
                                    2:getAFBC_v}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v}}

class Steady_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Steady_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Steady_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v()}
