from pyadh import *
from pyadh.default_p import *
from pyadh.ctransportCoefficients import smoothedHeaviside
from twp_step3d import *

coefficients = VOFCoefficients(LS_model=1,V_model=0,RD_model=3,ME_model=2,epsFact=epsFact_vof)

class Flat_H:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        signedDistance = (x[2] - self.waterLevel)
        return smoothedHeaviside(epsFact_consrv_heaviside*he,signedDistance)

H_init = Flat_H(waterLevel)
analyticalSolutions = None

def getDBC_vof(x,flag):
    if flag in outflow:
        return lambda x,t: 1.0
    elif flag == boundaryTags['upstream']:
        return H_init.uOfXT

dirichletConditions = {0:getDBC_vof}

initialConditions  = {0:H_init}

fluxBoundaryConditions = {0:'mixedFlow'}

def getAFBC_vof(x,flag):
    if flag in bottom+noflow:
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_vof}

diffusiveFluxBoundaryConditions = {0:{}}
