from pyadh import *
from pyadh.default_p import *
from pyadh.ctransportCoefficients import smoothedHeaviside
from threep_cylinder_md_2d import *

coefficients = VOFCoefficients(LS_model=2,
                               V_model=1,
                               RD_model=4,
                               ME_model=3,
                               epsFact=epsFact_vof)

class PerturbedSurface_H:
    def __init__(self,waterLevel,slopeAngle):
        self.waterLevel=waterLevel
        self.slopeAngle=slopeAngle
    def uOfXT(self,x,t):
        if useShock:
            return smoothedHeaviside(epsFact_consrv_heaviside*he,shockSignedDistance(x))
        else:
            surfaceNormal = [-sin(self.slopeAngle),cos(self.slopeAngle)]
            signedDistance = (x[0] - 0.5)*surfaceNormal[0]+(x[1] - self.waterLevel)*surfaceNormal[1]
            return smoothedHeaviside(epsFact_consrv_heaviside*he,signedDistance)

analyticalSolutions = None

def getDBC_vof(x,flag):
    if flag == domain.boundaryTags['top']:
        return lambda x,t: 1.0
    if flag == domain.boundaryTags['upstream']:
        #return lambda x,t: smoothedHeaviside(epsFact_consrv_heaviside*he,x[1]-waterLevel + 0.05*sin(pi*t))
        return lambda x,t: smoothedHeaviside(epsFact_consrv_heaviside*he,x[1]-waterLevel)
    if flag == domain.boundaryTags['downstream']:
        return lambda x,t: smoothedHeaviside(epsFact_consrv_heaviside*he,x[1]-waterLevel)

dirichletConditions = {0:getDBC_vof}

initialConditions  = {0:PerturbedSurface_H(waterLevel,0.0)}

#fluxBoundaryConditions = {0:'noFlow'}
fluxBoundaryConditions = {0:'mixedFlow'}

def getAFBC_vof(x,flag):
    if flag == boundaryTags['bottom']:
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_vof}

diffusiveFluxBoundaryConditions = {0:{}}
