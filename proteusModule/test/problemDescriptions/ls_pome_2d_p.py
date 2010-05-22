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

#coefficients = NCLevelSetCoefficients(V_model=4,RD_model=2,ME_model=5)
coefficients = NCLevelSetCoefficients(V_model=2,RD_model=0,ME_model=3)
#coefficients = CLevelSetCoefficients(V_model=4,RD_model=2,ME_model=5)
#coefficients = CLevelSetCoefficients(V_model=1,RD_model=-1,ME_model=2)
#coefficients = CLevelSetCoefficients(V_model=4,RD_model=-1,ME_model=5)

a=1.5
b=2.0
a=1.0
#b=1.5
b=1.0
class Bubble_phi:
    def __init__(self,center=[0.5,0.2],radius=1.0/8.0):
        self.radius  = radius
        self.center  = center
    def uOfX(self,X):
        dwaterLevel = X[1] - waterLevel
        return dwaterLevel
    #end
    def uOfXT(self,X,t):
        return self.uOfX(X)
    #end
#end Bubble_phi

class Tube_phi:
    def __init__(self,height):
        self.height=height
    def uOfX(self,X):
        xc = 0.5*L[0]
        yc = 0.75*L[1]
        if X[1] <= yc:
            return -(sqrt((X[1]-yc)**2 + (X[0]-xc)**2) - 0.5*L[0])
        else:
            return (X[1]-yc) -(sqrt((X[0]-xc)**2) - 0.5*L[0])
#         return X[1]-self.height
    def uOfXT(self,X,t):
        return self.uOfX(X)

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

dirichletConditions = {0:getDBC}
def getAFBC(x):
     pass

diffusiveFluxBoundaryConditions = {0:{}}
