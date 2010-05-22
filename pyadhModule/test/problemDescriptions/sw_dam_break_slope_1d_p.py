from pyadh import *
from pyadh.default_p import *
"""
 Dam Break over sloping bottom
"""
nd=1


L=(10.0,)
domain = Domain.RectangularDomain(L=L)
slope = 0.1
import math
def bathymetry(x):
    return slope*(L[0]-x[0])

def bathymetryGrad(x):
    return -slope
g =9.8
HR=0.0
HL=1.0

coefficients = ShallowWater(g=g,
                            nd=nd,
                            bathymetryFunc=bathymetry,
                            bathymetryGradientFunc=bathymetryGrad)

class DamBreakIC:
    def __init__(self,xc,HL,HR):
        self.xc = xc
        self.HL=HL
        self.HR=HR
    def uOfXT(self,x,t):
        if x[0] < self.xc:
            h = self.HL
        else:
            h = self.HR
        return h
class ZeroIC:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:DamBreakIC(0.5,HL,HR),
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


T=50.0


