from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from math import *
from clsvof import *

LevelModelType = CLSVOF.LevelModel
logEvent = Profiling.logEvent
name=soname

nd=2

coefficients = MyCoefficients(
    checkMass=checkMass,
    useMetrics=ct.useMetrics,
    timeOrder=ct.timeOrder,
    MONOLITHIC_CLSVOF=True,
    epsFactHeaviside=1.5,
    epsFactDirac=1.5,
    epsFactDiffusion=ct.lambdaFact,
    outputQuantDOFs=False,
    computeMetrics=True)
coefficients.variableNames=['u']

##################
# VELOCITY FIELD #
##################
def velx(X,t):
    x=X[0]
    y=X[1]
    #return -2*pi*(X[1]-0.5)
    return -2*np.sin(pi*x)**2*np.sin(pi*y)*np.cos(pi*y)*np.sin(2*pi*t/ct.T)

def vely(X,t):
    x=X[0]
    y=X[1]
    #return 2*pi*(X[0]-0.5)
    return 2*np.sin(pi*y)**2*np.sin(pi*x)*np.cos(pi*x)*np.sin(2*pi*t/ct.T)

velocityFieldAsFunction={0:velx, 
                         1:vely}

def exact_solution(X,t):
    radius = 0.15
    xc = 0.5
    yc = 0.75
    r = np.sqrt((X[0]-xc)**2 + (X[1]-yc)**2)
    return radius - r
exactSolution={0:exact_solution}
    
#####################
# INITIAL CONDITION #
#####################
class init_cond:
    def __init__(self,L):
        self.radius = 0.15
        self.xc=0.5
        self.yc=0.75
    def uOfXT(self,x,t):
        r = math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2)
        return self.radius - r
    
initialConditions  = {0:init_cond(L)}

#######################
# BOUNDARY CONDITIONS #
#######################
def getDBC(x,flag):
    None
dirichletConditions = {0:getDBC}

#Periodic Boundary Conditions
eps=1.0e-8
def getPDBC(x,flag):
    pass
periodicDirichletConditions = {0:getPDBC}

def zeroadv(x,flag):
    None
advectiveFluxBoundaryConditions =  {0:zeroadv}

fluxBoundaryConditions = {0:'outFlow'}
diffusiveFluxBoundaryConditions = {0:{}}
