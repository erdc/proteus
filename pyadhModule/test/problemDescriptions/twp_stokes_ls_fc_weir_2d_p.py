from pyadh import *
from pyadh.default_p import *

nd = 2
                     
analyticalSolution = None

polyfile = "weir"

coefficients = TwophaseStokes_LS(g=[0.0,9.8],nd=nd)
coefficients.eps=0.00001
coefficients.rho_2 = coefficients.rho_1
coefficients.nu_2 = coefficients.nu_1

flume_top = 0.61
flume_end = 21.4
inflow = 0.1
waterLevel = 0.5
jetHeight = waterLevel
weirHeight = .235
weirStart = 8.96
weirEnd = 8.98

def getDBC_phi_weir(x):
    #inflow
    if x[0] == flume_end:
        return lambda x,t: x[1]-waterLevel
#     #bottom
#     if x[1] == 0.0:
#         return lambda x,t: x[1]-waterLevel
#     #weir
#     if (x[0] >= weirStart and
#         x[0] <= weirEnd):
#         if x[1] <= weirHeight:
#             return lambda x,t: x[1]-waterLevel
#     #top
#     if x[1] == flume_top:
#         return lambda x,t: x[1]-waterLevel

def getDBC_p_weir(x):
    #top
    if x[1] == flume_top:
        return lambda x,t: 0.0

def getDBC_u_weir(x):
    #bottom
    if (x[1] == 0.0):
        return lambda x,t: 0.0
    #weir
    if (x[0] >= weirStart and
        x[0] <= weirEnd):
        if x[1] <= weirHeight:
            return lambda x,t: 0.0
    #inflow
    if (x[0] == flume_end):
        if x[1] <= waterLevel and x[1] >= waterLevel-jetHeight:
            return lambda x,t: -inflow
        else:
            return lambda x,t: 0.0

def getDBC_v_weir(x):
    #botom
    if (x[1] == 0.0):
        return lambda x,t: 0.0
    #weir
    if (x[0] >= weirStart and
        x[0] <= weirEnd):
        if x[1] <= weirHeight:
            return lambda x,t: 0.0
    #inflow
    if(x[0] == flume_end):
        return lambda  x,t: 0.0

dirichletConditions = {0:getDBC_phi_weir,
                       1:getDBC_p_weir,
                       2:getDBC_u_weir,
                       3:getDBC_v_weir}

def getAFBC_p_weir(x):
    #bottom
    if (x[1] == 0.0):
        return lambda x,t: 0.0
    #weir
    if (x[0] >= weirStart and
        x[0] <= weirEnd):
        if x[1] <= weirHeight:
            return lambda x,t: 0.0
    #inflow
    elif (x[0] == flume_end):
        if (x[1] <= waterLevel and x[1] >= waterLevel-jetHeight):
            return lambda x,t: -inflow
        else:
            return lambda x,t: 0.0

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow',
                          3:'outFlow'}

advectiveFluxBoundaryConditions =  {1:getAFBC_p_weir}

diffusiveFluxBoundaryConditions = {1:{}}

class SteadyNoWeir_phi:
    def __init__(self,waterLevel,inflow):
        self.waterLevel=waterLevel
        self.inflow=inflow
    def uOfXT(self,x,t):
        return x[1] - waterLevel

class SteadyNoWeir_p:
    def __init__(self,waterLevel,inflow):
        self.waterLevel=waterLevel
        self.inflow=inflow
    def uOfXT(self,x,t):
        if x[1] <= self.waterLevel:
            return ((flume_top-self.waterLevel)*coefficients.rho_2+(self.waterLevel-x[1])*coefficients.rho_1)*coefficients.g[1]
        else:
            return (flume_top-x[1])*coefficients.rho_2*coefficients.g[1]

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

initialConditions = {0:SteadyNoWeir_phi(waterLevel,inflow),
                     1:SteadyNoWeir_p(waterLevel,inflow),
                     2:SteadyNoWeir_u(waterLevel,inflow),
                     3:SteadyNoWeir_v(waterLevel,inflow)}

T=20.0
