from pyadh import *
from pyadh.default_p import *
from bubble import *
from math import *

coefficients = LevelSetCurvatureCoefficients(epsFact=epsFact_curvature,LSModel_index=2,nd=2)

fluxBoundaryConditions = {0:'outFlow'}

def getDBC(x,flag):
    pass

dirichletConditions = {0:getDBC}
a=1.5
b=2.0
a=1.0
#b=1.5
b=1.0
class Bubble_kappa:
    def __init__(self,center=[0.5,0.2],radius=1.0/8.0):
        self.radius  = radius
        self.center  = center
    def uOfX(self,X):
        dx = a*(X[0]-self.center[0]); dy = b*(X[1]-self.center[1])
        r = sqrt(dx**2 + dy**2)
        dBubble = self.radius - r
        dwaterLevel = X[1] - waterLevel
        if math.fabs(dBubble) < math.fabs(dwaterLevel):
            return 1.0/max(1.0e-8,r)
        else:
            return 0.0
        #return sqrt(dx**2 + dy**2) - self.radius
    #end
    def uOfXT(self,X,t):
        return self.uOfX(X)
    #end
#end Bubble_kappa

initialConditions  = {0:Bubble_kappa(center=[0.5*L[0],1.5*bubbleRadius],radius=bubbleRadius)}

def getAFBC(x,flag):
     pass
advectiveFluxBoundaryConditions =  {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{}}
