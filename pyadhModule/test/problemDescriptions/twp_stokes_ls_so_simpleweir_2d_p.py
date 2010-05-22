from pyadh import *
from pyadh.default_p import *

nd = 2
                     
analyticalSolution = None

polyfile = "simpleWeir"

coefficients = TwophaseStokes_LS_SO(g=[0.0,9.8],nd=nd,steady=False)
coefficients = TwophaseStokes_LS_SO(g=[0.0,9.8],nd=nd,steady=True)
coefficients.eps=1.0e-1
#coefficients.rho_1 = coefficients.rho_0
#coefficients.mu_1 = coefficients.mu_0

inflow = 0.1
waterLevel = 0.5
weirHeight = 0.4
weirStart = 1.4
weirEnd = 1.6
flumeEnd = 2.0
flumeTop = 1.0

def getDBC_p_weir(x):
    if x[1] == flumeTop:
        return lambda x,t: 0.0
        
def getDBC_u_weir(x):
    #bottom
    if (x[1] == 0.0):
        return lambda x,t: 0.0
    #weir
    if (x[0] >= weirStart and
        x[0] <= weirEnd):
        if (x[1] <= weirHeight):
            return lambda x,t: 0.0
    #inflow
    if (x[0] == 0.0):
        return lambda x,t: inflow

def getDBC_v_weir(x):
    #bottom
    if (x[1] == 0.0):
        return lambda x,t: 0.0
    #weir
    if (x[0] >= weirStart and
        x[0] <= weirEnd):
        if (x[1] <= weirHeight):
            return lambda x,t: 0.0
    #inflow
    if (x[0] == 0.0):
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p_weir,
                       1:getDBC_u_weir,
                       2:getDBC_v_weir}

def getAFBC_p_weir(x):
    #bottom
    if (x[1] == 0.0):
        return lambda x,t: 0.0
    #weir
    if (x[0] >= weirStart and
        x[0] <= weirEnd):
        if (x[1] <= weirHeight):
            return lambda x,t: 0.0
    #inflow
    elif (x[0] == 0.0):
        return lambda x,t: -inflow

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_weir}

diffusiveFluxBoundaryConditions = {0:{}}

class SteadyNoWeir_p:
    def __init__(self,waterLevel,inflow):
        self.waterLevel=waterLevel
        self.inflow=inflow
    def uOfXT(self,x,t):
        if x[1] > waterLevel:
            return (1.0-x[1])*coefficients.rho_1*coefficients.g[1]
        else:
            return ((1.0-waterLevel)*coefficients.rho_1 +
                    (waterLevel-x[1])*coefficients.rho_0
                    )*coefficients.g[1]

class SteadyNoWeir_u:
    def __init__(self,waterLevel,inflow):
        self.waterLevel=waterLevel
        self.inflow=inflow
    def uOfXT(self,x,t):
        return inflow

class SteadyNoWeir_v:
    def __init__(self,waterLevel,inflow):
        self.waterLevel=waterLevel
        self.inflow=inflow
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:SteadyNoWeir_p(waterLevel,inflow),
                     1:SteadyNoWeir_u(waterLevel,inflow),
                     2:SteadyNoWeir_v(waterLevel,inflow)}

T=10000.0
