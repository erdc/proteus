from pyadh import *
from pyadh.default_p import *
from math import *
from pyadh import RDLS
LevelModelType = RDLS.OneLevelRDLS
nd = 1


class ConstantVelocityCircleSquared:
    def __init__(self,radius=0.1,b=[1.],startX=0.25):
        self.radius = radius
        self.b      = b
        self.startX = startX
    def uOfXT(self,x,t):
        centerX = self.b[0]*t + self.startX
        return (x[0]-centerX)**2 - self.radius**2

b0 = numpy.array([1.],'d')
x0 = 0.5
r0 = 1./8.
analyticalSolution = {0:ConstantVelocityCircleSquared(r0,b0,x0)}

def Heaviside(p,eps=1.0e-2):
    if p < -eps:
        return 0.0
    if p > eps:
        return 1.0
    return 0.5*(1.+p/eps + 1./math.pi*math.sin(math.pi*p/eps))
def Se(p,eps=1.0e-2):
    return 2.0*(Heaviside(p,eps)-0.5)


coefficients = RedistanceLevelSet(epsFact=0.25,u0=ConstantVelocityCircleSquared(r0,b0,x0))


coefficients.variableNames=['u']
#zero level set?
if LevelModelType == RDLS.OneLevelRDLS:
    weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCs}
else:
    weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}


#now define the Dirichlet boundary conditions

def getDBC(x,flag):
    pass
    #if (x[1] == 0.0):
    #    return lambda x,t: 0.0
    #if (x[0] == 0.0 or
    #    x[0] == 1.0 or
    #    x[1] == 0.0 or
    #    x[1] == 1.0):
    #    return lambda x,t: 0.0
    
dirichletConditions = {0:getDBC}

initialConditions  = {0:analyticalSolution[0]}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 5.0e-1
