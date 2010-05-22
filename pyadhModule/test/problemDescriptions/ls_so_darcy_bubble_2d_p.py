from pyadh import *
from pyadh.default_p import *
from twp_darcy_ls_so_bubble_2d_p import nd

"""
The non-conservative level set description of a bubble in a two-phase Darcy flow
"""

##\ingroup test
#\file ls_so_darcy_bubble_2d_p.py
#
#\todo finish ls_so_darcy_bubble_2d_p.py doc

class DarcyNCLevelSetCoefficients(NCLevelSetCoefficients):
    def __init__(self,V_model=0,RD_model=-1,ME_model=1,EikonalSolverFlag=0):
        NCLevelSetCoefficients.__init__(self,V_model=V_model,RD_model=RD_model,ME_model=ME_model,
                                        EikonalSolverFlag=EikonalSolverFlag)
#end DarcyNC

coefficients = DarcyNCLevelSetCoefficients(V_model=0,RD_model=-1,ME_model=1,EikonalSolverFlag=0)


class Bubble_phi:
    def __init__(self,center=[0.5,0.2],radius=1.0/8.0):
        self.radius  =1.0/8.0
        self.center  =center
    def uOfX(self,X):
        dx = X[0]-self.center[0]; dy = X[1]-self.center[1]
        dist2 = dx**2 + dy**2 - self.radius**2
        #return dist2
        dist  = sqrt(dx**2 + dy**2) - self.radius
        return dist
        #if dist2 >= 0:
        #    return sqrt(dist2)
        #else:
        #    return -1.0*sqrt(-dist2)
        
    #end
    def uOfXT(self,X,t):
        return self.uOfX(X)
    #end
#end Bubble_phi
class Bubble_phiM:
    def __init__(self,center=[0.5,0.2],radius=1.0/8.0):
        self.radius  =1.0/8.0
        self.center  =center
    def uOfX(self,X):
        dx = X[0]-self.center[0]; dy = X[1]-self.center[1]
        dist2 = dx**2 + dy**2 - self.radius**2
        #return dist2
        dist  = sqrt(dx**2 + dy**2) - self.radius
        return -dist
        #if dist2 >= 0:
        #    return sqrt(dist2)
        #else:
        #    return -1.0*sqrt(-dist2)
        
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
#initialConditions  = {0:Bubble_phi(center=[0.5,0.2],radius=1.0/8.0)}
#bubble fall
initialConditions  = {0:Bubble_phiM(center=[0.5,0.8],radius=1.0/8.0)}


fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}



T = 4.0e-1
