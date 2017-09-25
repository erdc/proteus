from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from math import *
from oneD_advection import *

LevelModelType = NCLS.LevelModel
logEvent = Profiling.logEvent
name=soname+"_ls"

nd=2

def velx(X,t):
    return 1.0

def vely(X,t):
    return 0.0

velocityField={0:velx, 
               1:vely}

class init_cond:
    def __init__(self,L):
        self.radius = 0.15
        self.xc=0.5
        self.yc=0.75
    def uOfXT(self,x,t):
        import numpy as np
        beta = epsCoupez*he
        if problem==0:
            return -beta*np.tanh((x[0]-0.7)/beta)*np.tanh((x[0]-0.3)/beta)
        elif problem==1:
            if (x[0] >= 0.3 and x[0] <= 0.7):
                return 1.0
            else:
                return 0.0
        elif problem==2:
            if (x[0] <= 0.5):
                return (x[0]-0.3) 
            else:
                return -(x[0]-0.7)             
        else: #sin function
            return np.sin(2*np.pi*x[0])
            
RD_model=None
coefficients = MyCoefficients(
    epsCoupez=epsCoupez*he, 
    epsFact=epsFactHeaviside,
    checkMass=checkMass,
    RD_model=RD_model,
    useMetrics=useMetrics, 
    LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX, 
    STABILIZATION_TYPE=STABILIZATION_TYPE, 
    ENTROPY_TYPE=ENTROPY_TYPE, 
    cE=cE)

coefficients.variableNames=['u']
initialConditions  = {0:init_cond(L)}

#now define the Dirichlet boundary conditions
def getDBC(x,flag):
#    return lambda x,t: 0.
    None
dirichletConditions = {0:getDBC}

#Periodic Boundary Conditions
eps=1.0e-8
def getPDBC(x,flag):
    if x[0] <= eps or x[0] >= (L[0]-eps):
        return numpy.array([0.0,x[1],0.0])
periodicDirichletConditions = {0:getPDBC}
 
def zeroInflow(x):
    None
fluxBoundaryConditions = {0:'outFlow'}

def zeroadv(x,flag):
    return lambda x,t: 0.0

#advectiveFluxBoundaryConditions =  {}
advectiveFluxBoundaryConditions =  {0:zeroadv}
diffusiveFluxBoundaryConditions = {0:{}}
