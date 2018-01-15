from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from math import *
from clsvof import *

LevelModelType = CLSVOF.LevelModel
logEvent = Profiling.logEvent
name=soname

nd=3

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
    z=X[2]
    #return 1.0
    return 2*np.sin(pi*x)**2*np.sin(2*pi*y)*np.sin(2*pi*z)*np.cos(pi*t/3.0)

def vely(X,t):
    x=X[0]
    y=X[1]
    z=X[2]
#    return 1.0
    return -np.sin(2*pi*x)*np.sin(pi*y)**2*np.sin(2*pi*z)*np.cos(pi*t/3.0)

def velz(X,t):
    x=X[0]
    y=X[1]
    z=X[2]
#    return 1.0
    return -np.sin(2*pi*x)*np.sin(2*pi*y)*np.sin(pi*z)**2*np.cos(pi*t/3.0)

velocityFieldAsFunction={0:velx, 
                         1:vely,
                         2:velz}

def exact_solution(X,t):
    radius = 0.15
    xc = 0.35
    yc = 0.35
    zc = 0.35
    r = np.sqrt((X[0]-xc)**2 + (X[1]-yc)**2 + (X[2]-zc)**2)
    return radius - r
exactSolution={0:exact_solution}

#####################
# INITIAL CONDITION #
#####################
class init_cond:
    def __init__(self,L):
        self.radius = 0.15
        self.xc=0.35
        self.yc=0.35
        self.zc=0.35
        
    def uOfXT(self,X,t):
        r = math.sqrt((X[0]-self.xc)**2 + (X[1]-self.yc)**2 + (X[2]-self.zc)**2)
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
