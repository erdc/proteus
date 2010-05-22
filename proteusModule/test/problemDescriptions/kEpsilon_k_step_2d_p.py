from pyadh import *
from pyadh.default_p import *
from rans_step_2d import *
import step2d

"""
k-epsilon turbulence model for backward facing step
k equation
"""


class ConstantIC:
    def __init__(self,cval=0.0):
        self.cval=cval
    def uOfXT(self,x,t):
        return self.cval

coefficients = kEpsilon_k(flowModelID=0,
                          epsilonModelID=2,
                          nd=nd,
                          g=g,
                          nu=nu,
                          rho=rho)

kInflow = 0.003*inflow*inflow
#mwf hack
#coefficients.c_mu = 0.0#turn off nonlinearities

initialConditions = {0:ConstantIC(cval=kInflow*0.001)}

analyticalSolution = None


def getDBC_k(x,flag):
    if flag == boundaryTags['upstream']:
        return lambda x,t:kInflow
    
dirichletConditions = {0:getDBC_k}


fluxBoundaryConditions = {0:'outFlow'}

diffusiveFluxBoundaryConditions = {0:{}}
