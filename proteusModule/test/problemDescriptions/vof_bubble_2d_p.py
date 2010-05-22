from pyadh import *
from pyadh.default_p import *
from bubble import *

"""
The non-conservative level set description of a bubble in a two-phase flow
"""

##\ingroup test
#\file ls_bubble_2d_p.py
#
# \todo finish ls_bubble_2d_p.py

coefficients = VOFCoefficients(LS_model=2,V_model=1,RD_model=4,ME_model=3)

def smoothedHeaviside(eps,phi):
    import math
    if phi > eps:
        return 1.0
    elif phi <= - eps:
        return 0.0
    else:
        return 0.5*(1.0 + phi/eps + math.sin(math.pi*phi/eps)/math.pi);
class Bubble_H:
    def __init__(self,center=[0.5,0.2],radius=1.0/8.0):
        self.radius  = radius
        self.center  = center
    def uOfX(self,X):
        dx = X[0]-self.center[0]; dy = X[1]-self.center[1]
        dBubble = self.radius - sqrt(dx**2 + dy**2)
        dwaterLevel = X[1] - waterLevel
        if math.fabs(dBubble) < math.fabs(dwaterLevel):
            return smoothedHeaviside(epsFact_consrv_heaviside*L[0]/float(nnx-1),dBubble)
        else:
            return smoothedHeaviside(epsFact_consrv_heaviside*L[0]/float(nnx-1),dwaterLevel)
        #return sqrt(dx**2 + dy**2) - self.radius
    #end
    def uOfXT(self,X,t):
        return self.uOfX(X)
    #end
#end Bubble_H

analyticalSolutions = None

def getDBC(x,flag):
    pass

dirichletConditions = {0:getDBC}
#bubble rise
initialConditions  = {0:Bubble_H(center=[0.5*L[0],1.5*bubbleRadius],radius=bubbleRadius)}
initialConditions  = {0:Bubble_H(center=[bubbleCenter_x,bubbleCenter_y],radius=bubbleRadius)}
#initialConditions  = {0:Tube_phi(height=0.5*L[1])}
#bubble fall
#initialConditions  = {0:Bubble_H(center=[0.5,0.5],radius=bubbleRadius)}


#fluxBoundaryConditions = {0:'noFlow'}
fluxBoundaryConditions = {0:'outFlow'}

def getAFBC(x,flag):
    pass
advectiveFluxBoundaryConditions =  {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{}}
