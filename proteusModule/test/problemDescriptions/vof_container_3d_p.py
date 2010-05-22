from pyadh import *
from pyadh.default_p import *
from pyadh.ctransportCoefficients import smoothedHeaviside
from container import *
from pyadh import VOF

LevelModelType = VOF.OneLevelVOF
coefficients = VOFCoefficients(LS_model=1,
                               V_model=0,
                               RD_model=3,
                               ME_model=2,
                               epsFact=epsFact_vof)

class Flat_H:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        VOF=smoothedHeaviside(epsFact_consrv_heaviside*he,x[2]-waterLevel)
        return VOF

analyticalSolutions = None

def getDBC_vof(x,flag):
    if flag == domain.boundaryTags['left']:
        return lambda x,t: smoothedHeaviside(epsFact_consrv_heaviside*he,x[2]-waterLevel)
    if flag == domain.boundaryTags['right']:
        return lambda x,t: smoothedHeaviside(epsFact_consrv_heaviside*he,x[2]-waterLevel)

dirichletConditions = {0:getDBC_vof}

initialConditions  = {0:Flat_H(waterLevel)}

fluxBoundaryConditions = {0:'mixedFlow'}

def getAFBC_vof(x,flag):
    if flag in [boundaryTags['bottom'],boundaryTags['top'],boundaryTags['front'],boundaryTags['back'],boundaryTags['obstacle']]:
        return lambda x,t: 0.0
    if flag == 0:
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_vof}

diffusiveFluxBoundaryConditions = {0:{}}
