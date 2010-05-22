from pyadh import *
from pyadh.default_p import *
import math
nd = 2

L=(1.0,1.0,1.0)
#L=(math.pi,math.pi,1.0)
name = "ns_vortex_tsgs"
Re = 1.0e-6
nVortices=2
coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                             sigma=0.0,
                                             rho_0=1.0,nu_0=1.0/Re,
                                             rho_1=1.0,nu_1=1.0/Re,
                                             g=[0.0,0.0],
                                             nd=2,
                                             LS_model=None,
                                             KN_model=None,
                                             epsFact_density=None,
                                             stokes=False)

pSol = AnalyticalSolutions.VortexDecay_p(nVortices,Re)
uSol = AnalyticalSolutions.VortexDecay_u(nVortices,Re)
vSol = AnalyticalSolutions.VortexDecay_v(nVortices,Re)

analyticalSolution = {0:pSol,1:uSol,2:vSol}

def getDBC_pressure(x,flag):
    if (x[0] in [0.0,L[0]] or
        x[1] in [0.0,L[1]]):
        return lambda x,t: pSol.uOfXT(x,t)

def getDBC_u(x,flag):
    if (x[0] in [0.0,L[0]] or
        x[1] in [0.0,L[1]]):
        return lambda x,t: uSol.uOfXT(x,t)

def getDBC_v(x,flag):
    if (x[0] in [0.0,L[0]] or
        x[1] in [0.0,L[1]]):
        return lambda x,t: vSol.uOfXT(x,t)

# def isInflow(x):
#     if x[0] == 0.0 and uSol.uOfXT(x,0) >= 0:
#         return True
#     elif x[0] == L[0] and uSol.uOfXT(x,0) <= 0:
#         return True
#     elif x[1] == 0.0 and vSol.uOfXT(x,0) >= 0:
#         return True
#     elif x[1] == L[1] and vSol.uOfXT(x,0) <= 0:
#         return True
#     else:
#         return False

# def getDBC_pressure(x,flag):
#     if not isInflow(x):
#         return lambda x,t: pSol.uOfXT(x,t)

# def getDBC_u(x,flag):
#     if isInflow(x):
#         return lambda x,t: uSol.uOfXT(x,t)

# def getDBC_v(x,flag):
#     if isInflow(x):
#         return lambda x,t: vSol.uOfXT(x,t)

dirichletConditions = {0:getDBC_pressure,1:getDBC_u,2:getDBC_v}

initialConditions  = {0:pSol,1:uSol,2:vSol}

fluxBoundaryConditions = {0:'outFlow',1:'outFlow',2:'outFlow'}

# def getAFBC(x,flag):
#     pass
# def getDFBC(x,flag):
#     if not isInflow(x):
#         return lambda x,t: 0.0
# def getZeroBC(x,flag):
#     if not isInflow(x):
#         return lambda x,t: 0.0

# advectiveFluxBoundaryConditions =  {0:getAFBC,1:getAFBC,2:getAFBC}

# diffusiveFluxBoundaryConditions = {0:{},1:{1:getDFBC},2:{2:getDFBC}}

T = 0.01*Re

