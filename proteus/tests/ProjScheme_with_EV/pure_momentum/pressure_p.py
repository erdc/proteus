from math import *
from proteus import *
from proteus.default_p import *
from NS_convergence import *
from proteus.mprans import Pres

name = "pressure"

coefficients=Pres.Coefficients(modelIndex=PRESSURE_model,
                               fluidModelIndex=V_model,
                               pressureIncrementModelIndex=PINC_model,
                               useRotationalForm=False)

def getDBC_p(x,flag):
    return lambda x,t: 0.0

def getFlux(x,flag):
    return lambda x,t: 0.0

class getIBC_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:getIBC_p()}
dirichletConditions = {0:getDBC_p } # pressure bc are explicitly set
advectiveFluxBoundaryConditions = {0:getFlux}
