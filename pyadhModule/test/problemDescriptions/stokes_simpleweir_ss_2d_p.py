from pyadh import *
from pyadh.default_p import *

nd = 2

initialConditions = None

analyticalSolution = None

polyfile="simpleWeir"

coefficients = Stokes(g=[0.0,-9.8],nd=nd)

#now define the Dirichlet boundary conditions

flume_top = 2.0
flume_end = 20.0
inflow = 0.1
waterLevel = flume_top#single phase
weirHeight = 0.4
weirStart = 9.9
weirEnd = 10.1

def getDBC_p_weir(x):
    if x[1] == flume_top and x[0] == flume_end:
        return lambda x,t: 0.0

def getDBC_u_weir(x):
    if (x[1] == 0.0):
        return lambda x,t: 0.0
    if (x[0] >= weirStart and
        x[0] <= weirEnd):
        if x[1] <= weirHeight:
            return lambda x,t: 0.0
    if (x[0] == 0):
        return lambda x,t: inflow
def getDBC_v_weir(x):
    if (x[1] == 0.0):
        return lambda x,t: 0.0
    if (x[0] >= weirStart and
        x[0] <= weirEnd):
        if x[1] <= weirHeight:
            return lambda x,t: 0.0
    if(x[0] == 0):
        return lambda  x,t: 0.0
    #mwf what about a fixed lid?
    if (x[1] == flume_top):
        return lambda x,t:0.0
dirichletConditions = {0:getDBC_p_weir,
                       1:getDBC_u_weir,
                       2:getDBC_v_weir}

def getAFBC_p_weir(x):
    if (x[1] == 0.0):
        return lambda x,t: 0.0
    if (x[0] >= weirStart and
        x[0] <= weirEnd):
        if x[1] <= weirHeight:
            return lambda x,t: 0.0
    elif (x[0] == 0.0):
        return lambda x,t: -inflow
    if (x[1] == flume_top):
        return lambda x,t: 0.0

def getNoBC(x):
    pass
def getDFBC_u_weir(x):
    if (x[1] == flume_top):
        return 0.0 #try to get d/dn phi_t = 0
    if (x[0] == flume_end):
        return 0.0 #try to get d/dn phi_n = 0
    
fluxBoundaryConditions = {0:'outFlow',
                          1:'setFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_weir,1:getNoBC,2:getNoBC}

diffusiveFluxBoundaryConditions = {0:{},1:{1:getDFBC_u_weir},2:{}}

class SteadyNoWeir_p:
    def __init__(self,inflow):
        self.inflow=inflow
    def uOfXT(self,x,t):
        return (flume_top-x[1])*coefficients.rho*abs(coefficients.g[1])

class SteadyNoWeir_u:
    def __init__(self,inflow):
        self.inflow=inflow
    def uOfXT(self,x,t):
        return inflow

class SteadyNoWeir_v:
    def __init__(self,inflow):
        self.inflow=inflow
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:SteadyNoWeir_p(inflow),
                     1:SteadyNoWeir_u(inflow),
                     2:SteadyNoWeir_v(inflow)}

T=1.0
