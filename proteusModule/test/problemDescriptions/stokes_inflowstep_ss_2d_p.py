from pyadh import *
from pyadh.default_p import *

nd = 2
 
L=(1.0,1.0,1.0)
                    
analyticalSolution = None

#coefficients = Stokes(g=[0.0,-9.8],nd=nd)
#works with Dirichlet bc's on outflow (in press. eqn)
#and outflow on the momentum
coefficients = StokesP(g=[0.0,-9.8],nd=nd)
#works when set zero advective and diffusive flux on outflow (x=L)
#coefficients = StokesP(g=[0.0,0.0],nd=nd)

#now define the Dirichlet boundary conditions

inflow = 0.1
stepHeight = 0.5

def getDBC_p_weir(x):
    #mwf try to skip Dirichlet on p by using outflow condition
    pass
    #mwf orig, works ok
    #hydrostatic on outflow
    #if x[0] == L[0]:
    #    return lambda x,t: (1.0-x[1])*coefficients.rho*abs(coefficients.g[1])

def getDBC_u_weir(x):
    #bottom and top
    if (x[1] == 0.0 or x[1] == L[1]):
        return lambda x,t: 0.0
    #inflow/step
    if x[0] == 0.0:
        if x[1] <= stepHeight:
            return lambda x,t: 0.0
        else:
            return lambda x,t: inflow

def getDBC_v_weir(x):
    #bottom and top
    if (x[1] == 0.0 or x[1] == L[1]):
        return lambda x,t: 0.0
    #inflow/step
    if x[0] == 0.0:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p_weir,
                       1:getDBC_u_weir,
                       2:getDBC_v_weir}

def getAFBC_p_weir(x):
    #bottom and top
    if (x[1] == 0.0 or x[1] == L[1]):
        return lambda x,t: 0.0
    #inflow/step
    if x[0] == 0.0:
        if x[1] <= stepHeight:
            return lambda x,t: 0.0
        else:
            return lambda x,t: -inflow

def getAFBC_u_weir(x):
    if x[0] == L[0]:
        #if no gravity
        #return lambda x,t: 0.0
        #if gravity
        return lambda x,t: (1.0-x[1])*abs(coefficients.g[1])
def getDFBC_u_weir(x):
    if x[0] == L[0]:
        return lambda x,t: 0.0
    
#
#mwf orig also works with StokesP if set Dirichlet BC's in press. eqn
#fluxBoundaryConditions = {0:'outFlow',
#                          1:'outFlow',
#                          2:'outFlow'}
#works with StokesP if set advective and diffusive flux on outflow
fluxBoundaryConditions = {0:'outFlow',
                          1:'setFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_weir,1:getAFBC_u_weir}

diffusiveFluxBoundaryConditions = {0:{},1:{1:getDFBC_u_weir}}

class SteadyNoWeir_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return (1.0-x[1])*coefficients.rho*abs(coefficients.g[1])

class SteadyNoWeir_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return inflow

class SteadyNoWeir_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:SteadyNoWeir_p(),
                     1:SteadyNoWeir_u(),
                     2:SteadyNoWeir_v()}

T=20.0
