from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from math import *
from thelper_vof import *

LevelModelType = VOF.LevelModel
logEvent = Profiling.logEvent
name=soname

nd=2

coefficients = MyCoefficients(
    checkMass=checkMass,
    FCT=ct.FCT,
    LUMPED_MASS_MATRIX=ct.LUMPED_MASS_MATRIX, 
    STABILIZATION_TYPE=ct.STABILIZATION_TYPE, 
    ENTROPY_TYPE=ct.ENTROPY_TYPE, 
    cE=ct.cE, cK=ct.cK) 
coefficients.variableNames=['u']

##################
# VELOCITY FIELD #
##################
def velx(X,t):
    if ct.problem==0:
        return 1.0
    else:        
        return -2*pi*(X[1]-0.5)

def vely(X,t):
    if ct.problem==0:
        return 0.0
    else:
        return 2*pi*(X[0]-0.5)

velocityField={0:velx, 
               1:vely}

#####################
# INITIAL CONDITION #
#####################
class init_cond:
    def __init__(self,L):
        self.radius = 0.15
        self.xc=0.5
        self.yc=0.75
    def uOfXT(self,x,t):
        if ct.problem==0:
            if (x[0] >= 0.3 and x[0] <= 0.7):
                return 1.0
            else:
                return 0.0
        else: 
            r = math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2)
            slit = x[0] > self.xc-0.025 and x[0] < self.xc+0.025 and x[1] < self.yc+0.1125
            if (r<=self.radius and slit==False):
                return 1.
            else:
                return 0.
initialConditions  = {0:init_cond(L)}
analyticalSolution = {0:init_cond(L)}

#######################
# BOUNDARY CONDITIONS #
#######################
def getDBC(x,flag):
    if ct.problem==0:
        None
    else:
        return lambda x,t: 0.0
dirichletConditions = {0:getDBC}

#Periodic Boundary Conditions
eps=1.0e-8
def getPDBC(x,flag):
    if x[0] <= eps or x[0] >= (L[0]-eps):
        return numpy.array([0.0,x[1],0.0])
if ct.problem==0:
    periodicDirichletConditions = {0:getPDBC}

def zeroadv(x,flag):
    None
advectiveFluxBoundaryConditions =  {0:zeroadv}

fluxBoundaryConditions = {0:'outFlow'}
diffusiveFluxBoundaryConditions = {0:{}}
