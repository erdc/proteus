from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from math import *
from twoD_advection import *

LevelModelType = NCLS.LevelModel
logEvent = Profiling.logEvent
name=soname+"_ls"

nd=2

def velx(X,t):
    y = X[1]
    return -2*pi*(y-0.5)

def vely(X,t):
    x = X[0]
    return 2*pi*(x-0.5)

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
        scaling = 0.25
        if problem==0: # Rotation of dist function
            return scaling*(self.radius - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2))
        else: # Rotation of saturated dist function
            return 2*scaling*beta*(smoothedHeaviside(epsFactHeaviside*he,self.radius 
                                                     - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2))-0.5)

class exact_soln:
    def __init__(self,L):
        self.radius = 0.15
        self.xc=0.5
        self.yc=0.75
    def uOfXT(self,x,t):
        import numpy as np
        beta = epsCoupez*he
        if problem==0: # Rotation of dist function
            return self.radius - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2)
        else: # Rotation of saturated dist function
            return 2*beta*(smoothedHeaviside(epsFactHeaviside*he,self.radius 
                                             - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2))-0.5)

analyticalSolution = {0:exact_soln(L)}

RD_model=None
coefficients = MyCoefficients(
    epsCoupez=epsCoupez*he, 
    pure_redistancing=pure_redistancing,
    epsFactRedistancing=epsFactRedistance, #for signed function
    redistancing_tolerance=redist_tolerance,
    lambda_coupez=lambda_coupez,
    epsFact=epsFactHeaviside,
    checkMass=checkMass,
    RD_model=RD_model,
    useMetrics=useMetrics, 
    COUPEZ=COUPEZ,
    DO_REDISTANCING=DO_REDISTANCING, 
    SATURATED_LEVEL_SET=SATURATED_LEVEL_SET,
    STABILIZATION_TYPE=STABILIZATION_TYPE, 
    ENTROPY_TYPE=ENTROPY_TYPE)

coefficients.variableNames=['u']
initialConditions  = {0:init_cond(L)}

#now define the Dirichlet boundary conditions
def getDBC(x,flag):
    #return lambda x,t: 0.
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
