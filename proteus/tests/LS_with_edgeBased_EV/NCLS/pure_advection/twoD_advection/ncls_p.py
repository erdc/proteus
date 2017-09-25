from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from math import *
from twoD_advection import *
import numpy as np

LevelModelType = NCLS.LevelModel
logEvent = Profiling.logEvent
name=soname+"_ls"

nd=2

def velx(X,t):
    x = X[0]
    y = X[1]
    if (problem < 3): #rotation
        return -2*pi*(y-0.5)
    else:
        return -2*np.sin(pi*x)**2*np.sin(pi*y)*np.cos(pi*y)*np.sin(pi*t)

def vely(X,t):
    x = X[0]
    y = X[1]
    if (problem < 3): #rotation
        return 2*pi*(x-0.5)
    else: 
        return 2*np.sin(pi*y)**2*np.sin(pi*x)*np.cos(pi*x)*np.sin(pi*t)

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

        if problem==0: # Rotation of dist function
            return self.radius - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2)
        elif problem==1: # Rotation of smoothed disk
            return 2*beta*(smoothedHeaviside(epsFactHeaviside*he,self.radius - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2))-0.5)
        elif problem==2: # Rotation of Zalesak disk
            r = math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2)
            slit = x[0] > self.xc-0.025 and x[0] < self.xc+0.025 and x[1] < self.yc+0.1125
            if (r<=self.radius and slit==False):
                return beta
            else:
                return -beta
        else: # Periodic vortex
            r = math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2)
            if r<=self.radius:
                return beta
            else:
                return -beta
                
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
    None
dirichletConditions = {0:getDBC}

def zeroInflow(x):
    None

fluxBoundaryConditions = {0:'outFlow'}

def zeroadv(x,flag):
    return lambda x,t: 0.0

#advectiveFluxBoundaryConditions =  {}
advectiveFluxBoundaryConditions =  {0:zeroadv}
diffusiveFluxBoundaryConditions = {0:{}}
