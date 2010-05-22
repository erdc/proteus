from pyadh import *
from pyadh.default_p import *
nd=1


L=(10.0,1.0,1.0)
g = 9.8
H = 1.0
U = 10.0

coefficients = ShallowWater(g=g,
                            nd=nd)

class ConstantIC:
    def __init__(self,H):
        self.H=H
    def uOfXT(self,x,t):
        return self.H

class RarefactionsIC:
    def __init__(self,Lx,H,U):
        self.xc = Lx/2.0
        self.H=H
        self.U=U
    def uOfXT(self,x,t):
        import math
        if x[0]<self.xc:
            hu = -self.H*self.U
        else:
            hu = self.H*self.U
        return hu

initialConditions = {0:ConstantIC(H),
                     1:RarefactionsIC(L[0],H,U)}

def getDBC_h(x,flag):
    return None

def getDBC_hu(x,flag):
    if x[0] < 1.0e-8:
        return lambda x,t: -H*U
    elif x[0] > L[0] - 1.0e-8:
        return lambda x,t: H*U
    else:
        return None
    
dirichletConditions = {0:getDBC_h,
                       1:getDBC_hu}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow'}

def getAFBC_h(x,flag):
    return None
#     if (x[0] < 1.0e-8 or
#         x[0] > L[0] - 1.0e-8):
#         return lambda x,t: H*U
#     else:
#         return None
    
def getAFBC_hu(x,flag):
    return None

advectiveFluxBoundaryConditions =  {0:getAFBC_h,
                                    1:getAFBC_hu}

diffusiveFluxBoundaryConditions = {0:{},1:{}}

T=1.5


