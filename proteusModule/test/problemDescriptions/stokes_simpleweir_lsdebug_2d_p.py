from pyadh import *
from pyadh.default_p import *

nd = 2
                     
analyticalSolution = None

coefficients = Stokes(g=[0.0,9.8],nd=nd)

inflow = 0.1
waterLevel = 0.3
weirHeight = 0.1

def getDBC_phi_weir(x):
    #debug
    return lambda x,t: 0.0
#     #inflow
#     if x[0] == 0.0:
#         return lambda x,t: x[1]-0.5
#     #bottom
#     if x[1] == 0.0:
#         return lambda x,t: x[1]-0.5
#     #weir
#     if x[0] == 1.0:
#         if x[1] <= weirHeight:
#             return lambda x,t: x[1]-0.5
#     #top
#     if x[1] == 1.0:
#         return lambda x,t: 0.5

def getDBC_p_weir(x):
    #debug
    if x[1] == 0.0:
        if x[0] == 0.5:
            return lambda x,t: 0.0
#     #top
#     if x[1] == 1.0:
#         return lambda x,t: 0.0
        
def getDBC_u_weir(x):
    #debug
    if  x[1] == 1.0:
        if x[0] > 0.0 or x[1] < 1.0:
            return lambda x,t: 0.1
    if (x[0] == 0.0 or
        x[0] == 1.0 or
        x[1] == 0.0):
        return lambda x,t: 0.0
#     #bottom
#     if (x[1] == 0.0):
#         return lambda x,t: 0.0
#     #weir
#     if (x[0] == 1.0):
#         if (x[1] <= weirHeight):
#             return lambda x,t: 0.0
#     #inflow
#     if (x[0] == 0.0):
#         return lambda x,t: inflow

def getDBC_v_weir(x):
    if(x[1] == 1.0 or
       x[1] == 0.0 or
       x[0] == 1.0 or
       x[0] == 0.0):
        return lambda x,t: 0.0
#     #botom
#     if (x[1] == 0.0):
#         return lambda x,t: 0.0
#     #weir
#     if (x[0] == 1.0):
#         if (x[1] <= weirHeight):
#             return lambda x,t: 0.0
#     #inflow
#     if(x[0] == 0.0):
#         return lambda  x,t: 0.0

dirichletConditions = {0:getDBC_p_weir,
                       1:getDBC_u_weir,
                       2:getDBC_v_weir}

def getAFBC_p_weir(x):
    if x[0] == 0.0 or x[0] == 1.0 or x[1] == 0.0 or x[1] == 1.0:
        return 0.0
#     #bottom
#     if (x[1] == 0.0):
#         return lambda x,t: 0.0
#     #weir
#     elif (x[0] == 1.0):
#         if (x[1] <= weirHeight):
#             return lambda x,t: 0.0
#     #inflow
#     elif (x[0] == 0.0):
#         return lambda x,t: -inflow

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_weir}

diffusiveFluxBoundaryConditions = {0:{}}

# fluxBoundaryConditions = {0:'noFlow',
#                           1:'noFlow',
#                           2:'noFlow'}

# advectiveFluxBoundaryConditions =  {}

# diffusiveFluxBoundaryConditions = {0:{}}

class SteadyNoWeir_phi:
    def __init__(self,waterLevel,inflow):
        self.waterLevel=waterLevel
        self.inflow=inflow
    def uOfXT(self,x,t):
        return 0.0
        return x[1] - 0.5

class SteadyNoWeir_p:
    def __init__(self,waterLevel,inflow):
        self.waterLevel=waterLevel
        self.inflow=inflow
    def uOfXT(self,x,t):
        if x[1] <= self.waterLevel:
            return ((1.0-self.waterLevel)*coefficients.rho+(self.waterLevel-x[1])*coefficients.rho)*coefficients.g[1]
        else:
            return (1.0-x[1])*coefficients.rho*coefficients.g[1]

class SteadyNoWeir_u:
    def __init__(self,waterLevel,inflow):
        self.waterLevel=waterLevel
        self.inflow=inflow
    def uOfXT(self,x,t):
        return self.inflow

class SteadyNoWeir_v:
    def __init__(self,waterLevel,inflow):
        self.waterLevel=waterLevel
        self.inflow=inflow
    def uOfXT(self,x,t):
        return 0.0

# initialConditions = {0:SteadyNoWeir_phi(waterLevel,inflow),
#                      1:SteadyNoWeir_p(waterLevel,inflow),
#                      2:SteadyNoWeir_u(waterLevel,inflow),
#                      3:SteadyNoWeir_v(waterLevel,inflow)}
initialConditions = {0:SteadyNoWeir_p(waterLevel,inflow),
                     1:SteadyNoWeir_u(waterLevel,inflow),
                     2:SteadyNoWeir_v(waterLevel,inflow)}
initialConditions = None

#T=20.0
