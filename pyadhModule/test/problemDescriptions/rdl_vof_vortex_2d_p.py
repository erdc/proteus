from pyadh import *
from pyadh.default_p import *
from vortex import *
"""
The non-conservative level set description of a bubble in a two-phase flow
"""

##\ingroup test
#\file vof_vortex_2d_p.py
#
# \todo finish vof_vortex_2d_p.py

coefficients = VOFCoefficients(LS_model=0,V_model=0,RD_model=2,ME_model=1)

def Heaviside(phi):
    if phi > 0:
        return 1.0
    else:
        return 0.0
    
class Vortex_phi:
    def __init__(self,center=[0.5,0.75],radius=0.15):
        self.radius  = radius
        self.center  = center
    def uOfX(self,X):
        dx = X[0]-self.center[0]; dy = X[1]-self.center[1]
        dBubble = self.radius - sqrt(dx**2 + dy**2)
        return Heaviside(dBubble)
    #end
    def uOfXT(self,X,t):
        return self.uOfX(X)
    #end
#end Vortex_phi

analyticalSolutions = None

def getDBC(x):
    pass

dirichletConditions = {0:getDBC}

initialConditions  = {0:Vortex_phi()}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}
