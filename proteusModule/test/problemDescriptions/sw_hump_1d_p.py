from pyadh import *
from pyadh.default_p import *
nd=1


L=(10.0,1.0,1.0)
g = 9.8
H0=1.0
H0=0.0
HH=10.0
coefficients = ShallowWater(g=g,
                            nd=nd)

class HumpIC:
    def __init__(self,Lx,r,H,H0):
        self.xc = Lx/2.0
        self.r=r
        self.H=H
        self.H0=H0
    def uOfXT(self,x,t):
        import math
        if abs(x[0]-self.xc) < self.r:
            h = self.H*(1.0+cos(math.pi*(x[0]-self.xc)/self.r))+self.H0
        else:
            h = self.H0
        return h
class ZeroIC:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:HumpIC(L[0],0.1*L[0],HH,H0),
                      1:ZeroIC()}

def getDBC_h(x,flag):
    return None

def getDBC_u(x,flag):
    if x[0] in [0.0,L[0]]:
        return lambda x,t: 0.0
    else:
        return None

dirichletConditions = {0:getDBC_h,
                       1:getDBC_u}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow'}

def getAFBC_h(x,flag):
    if x[0] in [0.0,L[0]]:
        return lambda x,t: 0.0
    else:
        return None

def getAFBC_u(x,flag):
    return None

advectiveFluxBoundaryConditions =  {0:getAFBC_h,1:getAFBC_u}

diffusiveFluxBoundaryConditions = {0:{},1:{}}

T=10.0


