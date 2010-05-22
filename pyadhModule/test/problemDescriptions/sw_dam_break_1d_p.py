from pyadh import *
from pyadh.default_p import *
nd=1


L=(10.0,1.0,1.0)
g = 9.8
HR=0.0
HL=1.0
eddyViscosity = 0.0
coefficients = ShallowWater(g=g,
                            nd=nd,
                            eddyViscosity=eddyViscosity)

class DamBreakIC:
    def __init__(self,Lx,HL,HR):
        self.xc = Lx/2.0
        self.HL=HL
        self.HR=HR
    def uOfXT(self,x,t):
        import math
        if x[0]>self.xc:
            h = self.HL
        else:
            h = self.HR
        return h
class ZeroIC:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:DamBreakIC(L[0],HL,HR),
                      1:ZeroIC()}

def getDBC_h(x,flag):
    return None

def getDBC_hu(x,flag):
    if (x[0] < 1.0e-8 or
        x[0] > L[0] - 1.0e-8):
        return lambda x,t: 0.0
    else:
        return None
    
dirichletConditions = {0:getDBC_h,
                       1:getDBC_hu}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow'}

def getAFBC_h(x,flag):
    if (x[0] < 1.0e-8 or
        x[0] > L[0] - 1.0e-8):
        return lambda x,t: 0.0
    else:
        return None
    
def getAFBC_hu(x,flag):
    return None

advectiveFluxBoundaryConditions =  {0:getAFBC_h,
                                    1:getAFBC_hu}

diffusiveFluxBoundaryConditions = {0:{},1:{}}

T=20.0


