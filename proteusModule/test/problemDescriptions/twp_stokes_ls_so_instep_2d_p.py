from pyadh import *
from pyadh.default_p import *

nd = 2
 
L=(1.0,1.0,1.0)
                    
analyticalSolution = None

coefficients = TwophaseStokes_LS_SO(g=[0.0,9.8e-5],nd=nd)
coefficients.eps = 1.0e-8
coefficients.rho_1  = coefficients.rho_0
coefficients.nu_1 = coefficients.nu_1
#now define the Dirichlet boundary conditions

freeSurfaceHeight = 0.75
inflow = 0.1
stepHeight = 0.5

def getDBC_p_weir(x):
    #top
    if x[1] == L[1]:
        return lambda x,t: 0.0

def getDBC_u_weir(x):
    #bottom
    if x[1] == 0.0:
        return lambda x,t: 0.0
    #inflow/step
    if x[0] == 0.0:
        if x[1] <= stepHeight:
            return lambda x,t: 0.0
        else:
            return lambda x,t: inflow

def getDBC_v_weir(x):
    #bottom
    if x[1] == 0.0:
        return lambda x,t: 0.0
    #inflow/step
    if x[0] == 0.0:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p_weir,
                       1:getDBC_u_weir,
                       2:getDBC_v_weir}

def getAFBC_p_weir(x):
    #bottom
    if x[1] == 0.0:
        return lambda x,t: 0.0
    #inflow/step
    if x[0] == 0.0:
        if x[1] <= stepHeight:
            return lambda x,t: 0.0
        else:
            return lambda x,t: -inflow

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_weir}

diffusiveFluxBoundaryConditions = {0:{}}

class SteadyNoWeir_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if x[1] > freeSurfaceHeight:
            return (1.0-x[1])*coefficients.rho_1*coefficients.g[1]
        else:
            return ((1.0-freeSurfaceHeight)*coefficients.rho_1 +
                    (freeSurfaceHeight-x[1])*coefficients.rho_0
                    )*coefficients.g[1]

class SteadyNoWeir_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return -inflow

class SteadyNoWeir_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:SteadyNoWeir_p(),
                     1:SteadyNoWeir_u(),
                     2:SteadyNoWeir_v()}

T=20.0
