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

coefficients = NCLevelSetCoefficients(V_model=0,RD_model=-1,ME_model=1,EikonalSolverFlag=1)


class Bubble_phi:
    def __init__(self,center=[0.5,0.2],radius=1.0/8.0):
        self.radius  = radius
        self.center  = center
    def uOfX(self,X):
        dx = X[0]-self.center[0]; dy = X[1]-self.center[1]
        return self.radius - sqrt(dx**2 + dy**2)
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
initialConditions  = {0:Bubble_phi(center=[0.5*L[0],1.1*bubbleRadius],radius=bubbleRadius)}
#bubble fall
#initialConditions  = {0:Bubble_phi(center=[0.5,0.5],radius=bubbleRadius)}


fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
