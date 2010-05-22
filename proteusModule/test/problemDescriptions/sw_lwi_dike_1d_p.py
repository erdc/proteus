from pyadh import *
from pyadh.default_p import *

import lwi_dike_2dDomain
"""
 generate waves and runup on LWI dike from OPTICREST report


 units should be meters by default


 z=0         beach slope            back slope  
             (1/6)            B      (-1/3)
 |--------|----------------|-------|---------|
 0        x_bs            x_be    x_ce       x_bse=L

x_bs = inflow length (1 [m] default), 
z(x) = 0, [0,x_bs]

x_be = x_bs + beachLength (4.8 [m])
z(x) = 0.0 + (x-x_bs)beachSlope, [x_bs,x_be]

x_ce = x_be + B [m]
z(x) = z(x_be), [x_be,x_ce]

x_bse= x_ce + backSlopeLength
z(x) = z(x_ce) + (x-x_ce)*backSlope [x_ce,x_bse]

x_L  = x_bse [m]
z(x) = z(x_bse)

total domain height = z(x_bse) + domainHeightPad
"""

nd=1
#2d version of problem
lwi_domain = lwi_dike_2dDomain.lwi_dike_2d()

#describe bathymetry as above
L=(lwi_domain.x_L,)
bathymetry = lwi_domain.bathymetry
bathymetryGrad = lwi_domain.bathymetryGrad

testNumber = 2
if testNumber == 2:
    #still water depth
    h_swd = 0.7 #[m]
    #wave height
    H     = 0.117 #[m]
    #wave period
    T_s   = 2.446

else: #1 by default
    #still water depth
    h_swd = 0.7 #[m]
    #wave height
    H     = 0.155 #[m]
    #wave period
    T_s   = 1.959

domain = Domain.RectangularDomain(L=L)


g =9.82
HR=0.0
HL=1.0
bedFrictionCoefficient = 0.015
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
import math
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

T=7.*T_s


