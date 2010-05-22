from pyadh import *
from pyadh.default_p import *
"""
 generate waves and runup on an embankment




               
 z=0
             1/10 slope      1/s slope           1/s_2 slope
 |--------|----------------|-----------|-------|-----------|-----|
 0        x_bs            x_be         x_se   x_ce        x_bse  L

x_bs = 2 [m], 
z(x) = 0, [0,x_bs]

x_be = x_bs + 4.48 [m]
z(x) = 0.0 + (x-x_bs)1/10, [x_bs,x_be]

x_se = x_be + (h_c + h_s)*s, h_c and h_s from Sitanggang Lynnett report
z(x) = (x-x_be)1/s + z(x_be), [x_be,x_se] 

x_ce = x_se + B [m]
z(x) = z(x_se), [x_se,x_ce]

x_bse= x_ce + (h_c + h_s)*0.5*s_2
z(x) = z(x_ce) + (x-x_ce)*1/s_2 [x_ce,x_bse]

x_L  = x_bse + 0.5 [m]
z(x) = z(x_bse), [x_bse,L]

"""
nd=1

#describe bathymetry as above

x_bs = 2.0; x_be = x_bs + 4.48
h_c = 0.054; h_s = 0.081; s = 3.0 #Run 1 from Sitaggang report
x_se = x_be + (h_c+h_s)*s
B    = 0.0896;
x_ce = x_se + B
s_2  = -0.2  #back slope to avoid disc. bathymetry for now
x_bse= x_ce + (h_c+h_s)*0.5*abs(s_2)
x_L  = x_bse + 0.5
L=(x_L,)
#still water depth
h_swd = 0.529 #[m]
#wave height
H     = 0.107 #[m]
#wave period
T_s   = 1.549

domain = Domain.RectangularDomain(L=L)

import math
def bathymetry(x):
    if x[0] <= x_bs:  return 0.0
    if x[0] <= x_be:  return 0.0 + (x[0]-x_bs)*0.1 #beach slope
    if x[0] <= x_se:  return bathymetry([x_be]) + (x[0]-x_be)/s
    if x[0] <= x_ce:  return bathymetry([x_se])
    if x[0] <= x_bse: return bathymetry([x_ce]) + (x[0]-x_ce)/s_2
    return bathymetry([x_bse])
def bathymetryGrad(x):
    if x[0] <= x_bs:  return (0.0,)
    if x[0] <= x_be:  return (0.1,) #beach slope
    if x[0] <= x_se:  return (1.0/s,)
    if x[0] <= x_ce:  return (0.0,)
    if x[0] <= x_bse: return (1./s_2,)
    return (0.0,)

g =9.82
HR=0.0
HL=1.0
bedFrictionCoefficient = 0.02
bedFrictionPower       = 0.333
coefficients = ShallowWater(g=g,
                            nd=nd,
                            bedFrictionCoefficient=bedFrictionCoefficient,
                            bedFrictionPower=bedFrictionPower,
                            bathymetryFunc=bathymetry,
                            bathymetryGradientFunc=bathymetryGrad)

class ConstIC:
    def __init__(self,val=0.0):
        self.val=val
    def uOfXT(self,x,t):
        return self.val
class StillWaterIC:
    def __init__(self,val=0.0):
        self.val=val
    def uOfXT(self,x,t):
        return max(0.0,self.val-bathymetry(x))

initialConditions = {0:StillWaterIC(val=h_swd),
                     1:ConstIC(val=0.0)}

def getDBC_h(x,flag):
    if x[0] <= 1.0e-8:
        return lambda x,t:h_swd + H*math.sin(2.0*math.pi*t/T_s)

def getDBC_hu(x,flag):
    return None
    
dirichletConditions = {0:getDBC_h,
                       1:getDBC_hu}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow'}

def getAFBC_h(x,flag):
    return None
    
def getAFBC_hu(x,flag):
    return None

advectiveFluxBoundaryConditions =  {0:getAFBC_h,
                                    1:getAFBC_hu}

diffusiveFluxBoundaryConditions = {0:{},1:{}}

T=20.*T_s


