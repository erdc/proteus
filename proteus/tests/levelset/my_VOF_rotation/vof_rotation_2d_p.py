from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from rotation2D import *
from proteus.mprans import VOF
name=soname+"_vof"

"""
The non-conservative level set description of a bubble in a two-phase flow
"""

LevelModelType = VOF.LevelModel

coefficients = MyCoefficients(epsFact=epsFactHeaviside,checkMass=checkMass,useMetrics=useMetrics,ME_model=0,
                              EDGE_VISCOSITY=EDGE_VISCOSITY,
                              ENTROPY_VISCOSITY=ENTROPY_VISCOSITY,
                              FCT=FCT,
                              cK=cK,cE=cE,cMax=cMax)

def Heaviside(phi):
    if phi > 0:
        return 1.0
    elif phi < 0:
        return 0.0
    else:
        return 0.5

class init_cond:
    def __init__(self,center=[0.5,0.75,0.5],radius=0.15):
        self.radius  = radius
        self.center  = center
    def uOfX(self,X):
        dx = X[0]-self.center[0]; dy = X[1]-self.center[1];
        dBubble = self.radius - sqrt(dx**2 + dy**2)

        #Zalesak disk
        #dBubble = 1.0*(sqrt(dx**2 + dy**2) <= self.radius) - 1.0*(sqrt(dx**2 + dy**2) > self.radius)
        #xSlit1 = X[0] < self.center[0]+0.025
        #xSlit2 = X[0] > self.center[0]-0.025
        #xSlit = xSlit1 and xSlit2
        #ySlit = X[1] < 0.75+0.1125
        #slit = xSlit*ySlit
        #if (slit==1):
        #    dBubble = -1
        return smoothedHeaviside(epsFactHeaviside*he,dBubble)#Heaviside(dBubble)

        # SMOOTH 
        #sm = 0.25
        #dist2 = dx**2 + dy**2
        #return 0.5*(1-tanh(dist2/sm**2-1))

    def uOfXT(self,X,t):
        return self.uOfX(X)
    #end
#end init_cond

analyticalSolutions = None

def getDBC(x,flag):
    return lambda x,t: 1.0
    #pass

dirichletConditions = {0:getDBC}

initialConditions  = {0:init_cond(center=[0.0,0.5],radius=0.25)}
#initialConditions  = {0:init_cond(center=[0.5,0.75],radius=0.15)}

fluxBoundaryConditions = {0:'outFlow'}

#cek made no flux since v.n = 0 for this v
def getAFBC(x,flag):
   return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{}}
