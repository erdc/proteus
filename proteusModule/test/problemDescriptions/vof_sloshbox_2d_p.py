from pyadh import *
from pyadh.default_p import *
from pyadh.ctransportCoefficients import smoothedHeaviside
from sloshbox import *

"""
The non-conservative level set description of a bubble in a two-phase flow
"""

##\ingroup test
#\file ls_bubble_2d_p.py
#
# \todo finish ls_bubble_2d_p.py

coefficients = VOFCoefficients(LS_model=1,V_model=0,RD_model=3,ME_model=2,epsFact=epsFact_vof)

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
            #return (1.0-smoothedHeaviside(epsFact_consrv_heaviside*he,signedDistance))
            return smoothedHeaviside(epsFact_consrv_heaviside*he,signedDistance)

analyticalSolutions = None

def getDBC_vof(x,flag):
    if not closedTop:
        if x[1] == L[1]:
            return lambda x,t: 1.0
        else:
            return None
    else:
        return None

dirichletConditions = {0:getDBC_vof}

initialConditions  = {0:PerturbedSurface_H(waterLevel,slopeAngle)}

#fluxBoundaryConditions = {0:'noFlow'}
fluxBoundaryConditions = {0:'mixedFlow'}

def getAFBC_vof(x,flag):
    if closedTop:
        return lambda x,t: 0.0
    else:
        if x[1] == L[1]:
            return None
        else:
            return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC_vof}

diffusiveFluxBoundaryConditions = {0:{}}
