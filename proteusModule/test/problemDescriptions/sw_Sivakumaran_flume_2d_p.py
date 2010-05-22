from pyadh import *
from pyadh.default_p import *
"""
 Dam Break over smoothly varying topography taken from Charlies's thesis and Sivakumaran etal
"""
nd=2


L=(2.5,1.)
domain = Domain.RectangularDomain(L=L)

import math
def bathymetry(x):
    return 0.20*math.exp(-0.5*(((x[0]-L[0]*0.5)/0.24)**2))

def bathymetryGrad(x):
    return (-3.472*(x[0] - 1.25)*math.exp(-(8.6806*(x[0] - 1.25)**2)),0.0)
g =9.82
HR=0.0
HL=1.0
bedFrictionCoefficient = 0.016**2
bedFrictionPower       = 1.0/3.0
coefficients = ShallowWater(g=g,
                            nd=nd,
                            bedFrictionCoefficient=bedFrictionCoefficient,
                            bedFrictionPower=bedFrictionPower,
                            bathymetryFunc=bathymetry,
                            bathymetryGradientFunc=bathymetryGrad)
tailWaterHeight=0.1#0.1 0.2
class ConstIC:
    def __init__(self,val=0.0):
        self.val=val
    def uOfXT(self,x,t):
        return self.val

initialConditions = {0:ConstIC(val=tailWaterHeight),
                     1:ConstIC(val=0.0),
                     2:ConstIC(val=0.0)}

def getDBC_h(x,flag):
    if x[0] >= L[0]-1.0e-8:
        return lambda x,t:tailWaterHeight #m

def getDBC_hu(x,flag):
    if (x[0] < 1.0e-8):
        return lambda x,t: 0.03599
    else:
        return None
def getDBC_hv(x,flag):
    if (x[0] < 1.0e-8):
        return lambda x,t:0.0
    if (x[1] < 1.0e-8 or x[1] > L[1]-1.0e-8):
        return lambda x,t:0.0
    
dirichletConditions = {0:getDBC_h,
                       1:getDBC_hu,
                       2:getDBC_hv}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

def getAFBC_h(x,flag):
    if (x[0] < 1.0e-8):
        return lambda x,t: -0.03599
    
def getAFBC_hu(x,flag):
    return None
def getAFBC_hv(x,flag):
    return None

advectiveFluxBoundaryConditions =  {0:getAFBC_h,
                                    1:getAFBC_hu,
                                    2:getAFBC_hv}

diffusiveFluxBoundaryConditions = {0:{},1:{},2:{}}

T=15.0


