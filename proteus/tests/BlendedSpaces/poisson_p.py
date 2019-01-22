from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from math import *
from poisson import *

LevelModelType = BlendedSpaces.LevelModel
logEvent = Profiling.logEvent
name=soname

nd=2

coefficients = MyCoefficients(forceStrongConditions=False if PROBLEM_TYPE==2 else True, 
                              outputQuantDOFs=True,
                              epsilon=0.01,
                              he=1./(2.0*(nn-1)),
                              PROBLEM_TYPE=PROBLEM_TYPE,
                              AUTOMATED_ALPHA=AUTOMATED_ALPHA,
                              ALPHA_FOR_GALERKIN_SOLUTION=ALPHA_FOR_GALERKIN_SOLUTION)
coefficients.variableNames=['u']

#######################
# AUXILIARY FUNCTIONS #
#######################
def H(z,eps):
    return (z>=-eps)*(z<=eps)*np.maximum(0.0,0.5*(1+z/eps+1./pi*np.sin(pi*z/eps))) + 1.0*(z>eps)
    
def Hp(z,eps):
    return (z>=-eps)*(z<=eps)*0.5*(1+np.cos(pi*z/eps))/eps

def Hpp(z,eps):
    return (z>=-eps)*(z<=eps)*-pi/(2.*eps**2)*np.sin(pi*z/eps)

##################
# VELOCITY FIELD #
##################
def velx(X,t):
    return 1.0
def vely(X,t):
    return 3.0
velocityFieldAsFunction={0:velx,1:vely}

#####################
# BLENDING FUNCTION #
#####################
epsilon=(he/2.0)
def alpha(X,t):
    x=X[0]
    y=X[1]
    eps = 2*he    
    delta = 1.001*epsilon #1, 8, 14
    nB=4
    if PROBLEM_TYPE==0:
        return 0.5*(x==0.5) + 1.0*(x>0.5) # alpha=0,1 and left and right
    elif PROBLEM_TYPE==1:
        return 1.0
        return 1.0*(y>=eps)*(y<=0.8)
    else: #PROBLEM_TYPE==2
        alpha = (1.0-(1.0*(x>=4./9-0*delta)*(x<=5./9+0*delta)*(y>=4./9-0*delta)*(y<=5./9+0*delta))
                 -1*(y<=nB*delta) -1*(y>=1.-nB*delta)
                 -1*(x<=nB*delta)*(y>nB*delta)*(y<1.-nB*delta) -1*(x>=1.-nB*delta)*(y>nB*delta)*(y<1.-nB*delta))
        return alpha
        
blendingFunction={0:alpha}
alphaFunction={0:alpha}

###########################
# FORCE FIELD AS FUNCTION #
###########################
def force(X,t):
    x=X[0]
    y=X[1]
    if PROBLEM_TYPE==0:
        return 8*pi**2*np.sin(2*pi*x)*np.sin(2*pi*y)
    else: # PROBLEM_TYPE=1 or 2
        return 0.0
    
forceFieldAsFunction={0:force}

##################
# EXACT SOLUTION #
##################
# for 2D isotropic poisson eqn
def exact_soln(X,t):
    x=X[0]
    y=X[1]
    return np.sin(2*pi*x)*np.sin(2*pi*y)

def exact_gradX(X,t):
    x=X[0]
    y=X[1]
    return 2*pi*np.cos(2*pi*x)*np.sin(2*pi*y)

def exact_gradY(X,t):
    x=X[0]
    y=X[1]
    return 2*pi*np.sin(2*pi*x)*np.cos(2*pi*y)

class exact_solution:
    def uOfXT(self,X,t):
        return exact_soln(X,t)
    def duOfXT(self,X,t):
        return [exact_gradX(X,t),exact_gradY(X,t)]
    
analyticalSolution = {0:exact_solution()}
exactSolution={0:exact_soln}
exactGrad={0:exact_gradX,
           1:exact_gradY}

#####################
# INITIAL CONDITION #
#####################
# This is dummy since the problem is not time dependent
class init_cond:
    def __init__(self):
        pass
    def uOfXT(self,X,t):
        x=X[0]
        y=X[1]
        return 0.
initialConditions  = {0:init_cond()}

#######################
# BOUNDARY CONDITIONS #
#######################
eps=0.2*he

def getDBC(x,flag):
    if PROBLEM_TYPE==0:
        if (x[0]==0. or x[0]==1. or x[1]==0. or x[1]==1.):
            return lambda x,t: 0.0
    elif PROBLEM_TYPE==1:
        if x[1]==0:
            if x[0]<=1.0/3.0:
                return lambda x,t: 1.0
            else:
                return lambda x,t: 0.0        
        elif x[0]==0.:
            return lambda x,t: 1.0
        elif x[0]==1.:
            return lambda x,t: 0.0        
        elif x[1]==1:
            return lambda x,t: 0.0
    else: # this is really not used. I impose the BCs differently
        #    return lambda x,t: 1.
        if x[0]<=0+eps or x[0]>=1-eps or x[1]<=0+eps or x[1]>=1-eps:
            return lambda x,t: -1.0
        elif x[0] >= 4./9-eps and x[0] <= 5./9+eps and x[1] >= 4./9-eps and x[1] <= 5./9+eps: 
            return lambda x,t: 1.0
#           
dirichletConditions = {0:getDBC}

def zeroadv(x,flag):
    None
advectiveFluxBoundaryConditions =  {0:zeroadv}
fluxBoundaryConditions = {0:'outFlow'}
diffusiveFluxBoundaryConditions = {0:{}}
