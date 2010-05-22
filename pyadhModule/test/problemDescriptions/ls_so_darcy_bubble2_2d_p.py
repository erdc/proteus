from pyadh import *
from pyadh.default_p import *
from twp_darcy_ls_so_bubble_2d_p import nd

##\todo merge ls_so_darcy_bubble2 with ls_so_darcy_bubble

class DarcyNCLevelSetCoefficients(NCLevelSetCoefficients):
    def __init__(self):
        NCLevelSetCoefficients.__init__(self)
    def attachModels(self,modelList):
        self.q_v = modelList[self.flowModelIndex].q[('velocity',0)]
        self.ebq_v = modelList[self.flowModelIndex].ebq[('velocity',0)]
#end DarcyNC

coefficients = DarcyNCLevelSetCoefficients()


class Bubble_phi:
    def __init__(self,center=[0.5,0.2],radius=1.0/8.0):
        self.radius  =1.0/8.0
        self.center  =center
    def uOfX(self,X):
        dx = X[0]-self.center[0]; dy = X[1]-self.center[1]
        return sqrt(dx**2 + dy**2) - self.radius
        
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
initialConditions  = {0:Bubble_phi(center=[0.5,0.2],radius=1.0/8.0)}
#bubble fall
#initialConditions  = {0:Bubble_phi(center=[0.5,0.8],radius=1.0/8.0)}


fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}



T = 1.0e1
