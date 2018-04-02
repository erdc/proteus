from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from math import *
from clsvof import *

LevelModelType = CLSVOF.LevelModel
logEvent = Profiling.logEvent
name=soname

nd=2
if ct.test_case>2:
    nd=3

coefficients = MyCoefficients(
    useMetrics=useMetrics,
    timeOrder=timeOrder,
    epsFactHeaviside=1.5,
    epsFactDirac=1.5,
    lambdaFact=lambdaFact,
    outputQuantDOFs=False,
    computeMetrics=computeMetrics)
coefficients.variableNames=['u']

##################
# VELOCITY FIELD #
##################
def velx(X,t):
    x=X[0]
    y=X[1]
    if ct.test_case==1:
        return -2*np.sin(pi*x)**2*np.sin(pi*y)*np.cos(pi*y)*np.sin(2*pi*t/T)
    elif ct.test_case==2 or ct.test_case==3:
        return -2*pi*(y-0.5)
    else:
        z=X[2]
        return 2*np.sin(pi*x)**2*np.sin(2*pi*y)*np.sin(2*pi*z)*np.cos(pi*t/3.0)

def vely(X,t):
    x=X[0]
    y=X[1]
    if ct.test_case==1:
        return 2*np.sin(pi*y)**2*np.sin(pi*x)*np.cos(pi*x)*np.sin(2*pi*t/T)
    elif ct.test_case==2 or ct.test_case==3:
        return 2*pi*(x-0.5)
    else:
        z=X[2]
        return -np.sin(2*pi*x)*np.sin(pi*y)**2*np.sin(2*pi*z)*np.cos(pi*t/3.0)

def velz(X,t):
    x=X[0]
    y=X[1]
    z=X[2]
    if ct.test_case==3:
        return 0.
    else:
        return -np.sin(2*pi*x)*np.sin(2*pi*y)*np.sin(pi*z)**2*np.cos(pi*t/3.0)

if ct.test_case==1 or ct.test_case==2:
    velocityFieldAsFunction={0:velx,
                             1:vely}
else:
    velocityFieldAsFunction={0:velx,
                             1:vely,
                             2:velz}

##################
# EXACT SOLUTION #
##################
def exact_solution(X,t):
    radius = 0.15
    if nd==2:
        xc = 0.5
        yc = 0.75
        r = np.sqrt((X[0]-xc)**2 + (X[1]-yc)**2)
        return radius - r
    else: #3D
        if ct.test_case==3: #solid rotation
            xc = 0.5
            yc = 0.75
            zc = 0.25
        else: #LeVeque test
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
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return exact_solution(x,0)
initialConditions  = {0:init_cond()}

#######################
# BOUNDARY CONDITIONS #
#######################
def getDBC(x,flag):
    return lambda x,t: -1.0
dirichletConditions = {0:getDBC}

def zeroadv(x,flag):
    None
advectiveFluxBoundaryConditions =  {0:zeroadv}

fluxBoundaryConditions = {0:'outFlow'}
diffusiveFluxBoundaryConditions = {0:{}}
