from pyadh import *
from pyadh.default_p import *
from pome import *

"""
The non-conservative level set description of a bubble in a two-phase flow
"""

##\ingroup test
#\file ls_bubble_2d_p.py
#
# \todo finish ls_bubble_2d_p.py

coefficients = CLevelSetCoefficients(V_model=2,RD_model=0,ME_model=4)

waterLevel = 2.0*L[1]#0.9*L[1]

def Heaviside(phi):
    if phi > 0:
        return 1.0
    else:
        return 0.0
class Bubble_phi:
    def __init__(self,center=[0.5,0.2],radius=1.0/8.0):
        self.radius  = radius
        self.center  = center
    def uOfX(self,X):
        dx = X[0]-self.center[0]; dy = X[1]-self.center[1]
        dBubble = self.radius - sqrt(dx**2 + dy**2)
        dwaterLevel = X[1] - waterLevel
        if math.fabs(dBubble) < math.fabs(dwaterLevel):
            return Heaviside(dBubble)
        else:
            return Heaviside(dwaterLevel)
        #return sqrt(dx**2 + dy**2) - self.radius
    #end
    def uOfXT(self,X,t):
        return self.uOfX(X)
    #end
#end Bubble_phi

analyticalSolutions = None

def getDBC(x):
    pass

dirichletConditions = {0:getDBC}
#bubble rise
initialConditions  = {0:Bubble_phi(center=[0.5*L[0],1.5*bubbleRadius],radius=bubbleRadius)}
initialConditions  = {0:Bubble_phi(center=[bubbleCenter_x,bubbleCenter_y],radius=bubbleRadius)}
#initialConditions  = {0:Tube_phi(height=0.5*L[1])}
#bubble fall
#initialConditions  = {0:Bubble_phi(center=[0.5,0.5],radius=bubbleRadius)}


#fluxBoundaryConditions = {0:'noFlow'}
fluxBoundaryConditions = {0:'outFlow'}

def getAFBC(x):
    pass
advectiveFluxBoundaryConditions =  {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{}}
