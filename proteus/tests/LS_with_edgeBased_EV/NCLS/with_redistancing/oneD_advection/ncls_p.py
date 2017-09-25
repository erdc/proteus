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
        pass
    def uOfXT(self,x,t):
        import numpy as np
        beta = epsCoupez*he
        scaling = 0.25
        if problem == 0: # distance function
            if (x[0] <= 0.5):
                return scaling*(x[0]-0.3) 
            else:
                return -scaling*(x[0]-0.7) 
        else: # Saturated distance function
            return -beta*scaling*np.tanh((x[0]-0.7)/beta)*np.tanh((x[0]-0.3)/beta)
    
class exact_soln:
    def __init__(self,L):
        pass
    def uOfXT(self,x,t):
        import numpy as np
        beta = epsCoupez*he
        if problem == 0: # distance function
            if (x[0] <= 0.5):
                return (x[0]-0.3) 
            else:
                return -(x[0]-0.7) 
        else: # Saturated distance function
            return -beta*np.tanh((x[0]-0.7)/beta)*np.tanh((x[0]-0.3)/beta)

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
    DO_SMOOTHING=DO_SMOOTHING,
    DO_REDISTANCING=DO_REDISTANCING, 
    SATURATED_LEVEL_SET=SATURATED_LEVEL_SET,
    STABILIZATION_TYPE=STABILIZATION_TYPE, 
    ENTROPY_TYPE=ENTROPY_TYPE)

coefficients.variableNames=['u']
initialConditions  = {0:init_cond(L)}

#now define the Dirichlet boundary conditions
def getDBC(x,flag):
    None
dirichletConditions = {0:getDBC}

#Periodic Boundary Conditions
eps=1.0e-8
def getPDBC(x,flag):
    if x[0] <= eps or x[0] >= (L[0]-eps):
        return numpy.array([0.0,x[1],0.0])
periodicDirichletConditions = {0:getPDBC}
 
def zeroInflow(x):
    return lambda x,t: 0.0
fluxBoundaryConditions = {0:'outFlow'}

def zeroadv(x,flag):
    return lambda x,t: 0.0

#advectiveFluxBoundaryConditions =  {}
advectiveFluxBoundaryConditions =  {0:zeroadv}
diffusiveFluxBoundaryConditions = {0:{}}
