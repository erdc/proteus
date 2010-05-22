from pyadh import *
from pyadh.default_p import *
nd=2
L=(10.0,10.0,1.0)
g = 9.8
H = 1.0
H0= 0.5

eddyViscosity = 0.0
coefficients = ShallowWater(g=g,
                            nd=nd,
                            eddyViscosity=eddyViscosity)

class DamBreakIC:
    def __init__(self,Lx,Ly,M,HL,HR):
        self.xc = Lx/2.0
        self.yc = Ly/2.0
        self.M=M
        self.HL=HL
        self.HR=HR
    def uOfXT(self,x,t):
        import math
        if (x[0] - self.xc)*self.M + (x[1] - self.yc) < 0:
            h = self.HL
        else:
            h = self.HR
        return h
class ZeroIC:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:DamBreakIC(L[0],L[1],-0.5,H,H0),
                     1:ZeroIC(),
                     2:ZeroIC()}

def getDBC_h(x,flag):
    return None

def getDBC_u(x,flag):
    if x[0] in [0.0,L[0]]:
        return lambda x,t: 0.0
    else:
        return None

def getDBC_v(x,flag):
    if x[1] in [0.0,L[1]]:
        return lambda x,t: 0.0
    else:
        return None

dirichletConditions = {0:getDBC_h,
                       1:getDBC_u,
                       2:getDBC_v}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

def getAFBC_h(x,flag):
    if (x[0] in [0.0,L[0]] or
        x[1] in [0.0,L[1]]):
        return lambda x,t: 0.0
    else:
        return None

def getAFBC_u(x,flag):
    return None

def getAFBC_v(x,flag):
    return None

advectiveFluxBoundaryConditions =  {0:getAFBC_h,1:getAFBC_u,2:getAFBC_v}

diffusiveFluxBoundaryConditions = {0:{},1:{},2:{}}

T=1.0#10.0


