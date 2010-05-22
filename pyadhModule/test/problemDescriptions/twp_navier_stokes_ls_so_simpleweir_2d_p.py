from math import log,exp
from pyadh import *
from pyadh.default_p import *
from simpleweir import *
"""
Two-phase incompressible Navier-Stokes flow in a flume with a rectangular weir
"""

## \page Tests Test Problems
#\ref twp_navier_stokes_ls_so_simpleweir_2d_p.py "Two-phase incompressible Navier-Stokes flow in a flume with a rectangular weir"
# 

##\ingroup test
# \file twp_navier_stokes_ls_so_simpleweir_2d_p.py
# @{
#
# \brief Two-phase incompressible Navier-Stokes flow over a weir.
#
# The governing equations are described by the
# pyadh::TransportCoefficients::TwophaseNavierStokes_LS_SO class. The
# domain and initial/boundary conditions are slip conditions and no
# flow on the lower boundary, constant pressure/outflow on the top boundary,
# outflow on the downstream boundary, and prescribed velocity on the
# inflow boundary.
#\image html simpleweir_1.jpg
#\image html simpleweir_2.jpg
#
# <A href="https://adh.usace.army.mil/pyadh-images/simpleweir.avi"> AVI Animation of velocity and pressure </A>
#

nd = 2
                     
analyticalSolution = None


coefficients = TwophaseNavierStokes_LS_SO(epsFact=epsFact,
                                          rho_0 = rho_0,
                                          nu_0 = nu_0,
                                          rho_1 = rho_1,
                                          nu_1 = nu_1,
                                          g=g,
                                          nd=nd)


p01 = log(0.01)

def flowRate(t):
    #ramp up flow rate
    #when two flume volumes would have passed through at the inlow rate we're at 99%
    return inflow*(1.0 - exp( p01*(t / (2.0*flumeEnd/inflow))))
def getDBC_p_weir(x):
    #top closed
    #if x[0] == 0.0:
    #    if x[1] == L[1]:
    #        return lambda x,t: 0.0
    #top open
    if x[1] == flumeTop:
        return lambda x,t: 0.0
#     if x[1] == 0:
#         if x[0] >= waterLevel:
#             return lambda x,t: -((flumeTop - waterLevel)*rho_1 + waterLevel*rho_0)*g[1]

def getDBC_u_weir(x):
    #if x[1] >= (flumeTop-1.0e-8):
    #    return lambda x,t: 0.0
    #bottom
    #if x[1] == 0.0:
    #    return lambda x,t: 0.0
    #weir
    #if (x[0] >= weirStart and
    #    x[0] <= weirEnd):
    #    if (x[1] <= weirHeight):
    #        return lambda x,t: 0.0
    #inflow
#     if x[0] == 0.0:
#         if x[1] <= waterLevel:
#             return lambda x,t: inflow
#         else:
#             pass
#             #return lambda x,t: 0.0
#     if x[0] >= flumeEnd - 1.0e-8:
#         if x[1] <= waterLevel:
#             return lambda x,t: inflow
#         else:
#             pass
#             #return lambda x,t: 0.0
    if x[0] in [0.0,flumeEnd]:
        return lambda x,t: 0.0

def getDBC_v_weir(x):
    if x[1] == 0.0:
        if x[0] <= 1.0:
            return lambda x,t: flowRate(t)#inflow
        if x[0] >= 19.0:
            return lambda x,t: -flowRate(t)#inflow
    #bottom
    #if x[1] == 0.0:
    #    return lambda x,t: 0.0
    #top
    #if x[1] == flumeTop:
    #    return lambda x,t: 0.0
    #weir
    #if (x[0] >= weirStart and
    #    x[0] <= weirEnd):
    #    if (x[1] <= weirHeight):
    #        return lambda x,t: 0.0
    #inflow
#     if x[0] == 0.0:
#         return lambda x,t: 0.0
#     if x[0] >= flumeEnd - 1.0e-8:
#         return lambda x,t: 0.0

dirichletConditions = {0:getDBC_p_weir,
                       1:getDBC_u_weir,
                       2:getDBC_v_weir}

def getAFBC_p_weir(x):
    if x[1] == flumeTop:
        pass
    elif x[1] == 0.0:
        if x[0] <= 1.0:
            return lambda x,t: -flowRate(t)#inflow
        if x[0] >= 19.0:
            return lambda x,t: flowRate(t)#inflow
        else:
            return lambda x,t: 0.0
    else:
        return lambda x,t: 0.0
#     if (x[0] >= (flumeEnd - 1.0e-8) or
#         x[1] >= (flumeTop -1.0e-8)):
#         pass
#    elif x[0] == 0.0:
#        if x[1] <= waterLevel:
#            return lambda x,t: -inflow
#        else:
#            return lambda x,t: 0.0
#    elif x[0] >= (flumeEnd-1.0e-8):
#        if x[1] <= waterLevel:
#            return lambda x,t: inflow
#        else:
#            return lambda x,t: 0.0

def getAFBC_u_weir(x):
    pass

def getAFBC_v_weir(x):
    pass

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_weir,
                                    1:getAFBC_u_weir,
                                    2:getAFBC_v_weir}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{},
                                   2:{}}

class SteadyNoWeir_p:
    def __init__(self,waterLevel,inflow):
        self.waterLevel=waterLevel
        self.inflow=inflow
    def uOfXT(self,x,t):
        if x[1] > waterLevel:
            return -(flumeTop-x[1])*coefficients.rho_1*coefficients.g[1]
        else:
            return -((flumeTop-waterLevel)*coefficients.rho_1 +
                    (waterLevel-x[1])*coefficients.rho_0
                    )*coefficients.g[1]

class SteadyNoWeir_u:
    def __init__(self,waterLevel,inflow):
        self.waterLevel=waterLevel
        self.inflow=inflow
    def uOfXT(self,x,t):
        return 0.0
#         return inflow

class SteadyNoWeir_v:
    def __init__(self,waterLevel,inflow):
        self.waterLevel=waterLevel
        self.inflow=inflow
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:SteadyNoWeir_p(waterLevel,inflow),
                     1:SteadyNoWeir_u(waterLevel,inflow),
                     2:SteadyNoWeir_v(waterLevel,inflow)}

