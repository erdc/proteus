from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from math import *
from ncls import *

LevelModelType = NCLS.LevelModel
logEvent = Profiling.logEvent
name=soname

nd=2

RD_model=None
coefficients = MyCoefficients(
    checkMass=checkMass,
    RD_model=RD_model,
    epsFact=epsFactHeaviside, 
    # Choice of numerical method #
    STABILIZATION_TYPE=ct.STABILIZATION_TYPE, 
    COUPEZ=ct.COUPEZ,
    DO_REDISTANCING=ct.DO_REDISTANCING, 
    SATURATED_LEVEL_SET=ct.SATURATED_LEVEL_SET,
    ENTROPY_TYPE=ct.ENTROPY_TYPE, 
    # Coupez and re-distancing parameters #
    pure_redistancing=ct.pure_redistancing,
    epsCoupez=epsCoupez*he, 
    epsFactRedistancing=epsFactRedistance, #for signed function
    redistancing_tolerance=redist_tolerance,
    lambda_coupez=lambda_coupez,
    # Entropy viscosity parameter #
    cE=ct.cE)
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
    def __init__(self,L,scaling=0.25):
        self.radius=0.15
        self.xc=0.5
        self.yc=0.75
        self.scaling=scaling
    def uOfXT(self,x,t):
        import numpy as np
        beta = epsCoupez*he
        if ct.problem==0: #1D problem
            if ct.level_set_function == 0: # distance function
                if (x[0] <= 0.5):
                    return self.scaling*(x[0]-0.3) 
                else:
                    return -self.scaling*(x[0]-0.7) 
            else: # Saturated distance function
                return -beta*self.scaling*np.tanh((x[0]-0.7)/beta)*np.tanh((x[0]-0.3)/beta)
        else: #2D problem
            if ct.level_set_function == 0: #distance function 
                return self.scaling*(self.radius - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2))
            else: #Saturated distance function
                return 2*beta*self.scaling*(smoothedHeaviside(epsFactHeaviside*he,self.radius 
                                                              - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2))-0.5)
                
    
analyticalSolution = {0:init_cond(L,scaling=1.0)}
initialConditions  = {0:init_cond(L,scaling=0.25 if ct.COUPEZ or ct.DO_REDISTANCING or ct.pure_redistancing else 1.0)}

#now define the Dirichlet boundary conditions
def getDBC(x,flag):
    None
dirichletConditions = {0:getDBC}

#Periodic Boundary Conditions
eps=1.0e-8
def getPDBC(x,flag):
    if x[0] <= eps or x[0] >= (L[0]-eps):
        return numpy.array([0.0,x[1],0.0])
if ct.problem==0:
    periodicDirichletConditions = {0:getPDBC}
 
def zeroInflow(x):
    return lambda x,t: 0.0
fluxBoundaryConditions = {0:'outFlow'}

def zeroadv(x,flag):
    return lambda x,t: 0.0

#advectiveFluxBoundaryConditions =  {}
advectiveFluxBoundaryConditions =  {0:zeroadv}
diffusiveFluxBoundaryConditions = {0:{}}
