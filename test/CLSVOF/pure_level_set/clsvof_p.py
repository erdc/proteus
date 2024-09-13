from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from math import *
try:
    from .clsvof import *
except:
    from clsvof import *

LevelModelType = CLSVOF.LevelModel
logEvent = Profiling.logEvent
name=soname

nd=2
if ct.test_case>2:
    nd=3

epsFactHeaviside=1.5
if ct.test_case==4:
    epsFactHeaviside=0.5

coefficients = MyCoefficients(useMetrics=useMetrics,
                              doSpinUpStep=doSpinUpStep,
                              epsFactHeaviside=epsFactHeaviside,
                              epsFactDirac=epsFactHeaviside,
                              lambdaFact=lambdaFact,
                              outputQuantDOFs=True,
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
def zalesak_disk_per_point(X,t):
    xc = 0.5
    yc = 0.75
    radius=0.15
    x=X[0]
    y=X[1]
    # distance to center of the disk
    r = np.sqrt((x-xc)**2 + (y-yc)**2)
    dist = radius - r # distance to circle
    # coordinates to slot coorners
    xslot1 = xc-0.025
    xslot2 = xc+0.025
    yslot1 = 0.75 - np.sqrt(radius**2-0.025**2)
    yslot2 = 0.85

    #distance to the boundary of the slot
    aux1 = np.abs(x-xslot1)
    aux2 = np.abs(x-xslot2)
    aux3 = np.abs(y-yslot2)
    aux4 = np.abs(y-yslot1)

    if (y > yslot1): #above disk
        if (xslot1 < x and x <= xslot2 and y <= yslot2): #inside slot
            return -np.min([aux1,aux2,aux3])
        else: #Not inside slot
            if x <= xslot1: #left of slot
                if y <= yslot2: #down top boundary of slot
                    return np.min([dist,aux1])
                else: #above top boundary of slot
                    return np.min([dist,np.sqrt(aux1**2+aux3**2)])
            elif x >= xslot2: #right of slot
                if y <= yslot2: #down top boundary of slot
                    return np.min([dist,aux2])
                else: #above top boundary of slot
                    return np.min([dist,np.sqrt(aux2**2+aux3**2)])
            else: #x-coordiante is within slot
                return np.min([dist,aux3])
    else: #below disk
        if x > xslot1 and x < xslot2: #x-coordinate is within slot
            return  -np.min([np.sqrt(aux1**2 + aux4**2), np.sqrt(aux2**2 + aux4**2)])
        else:
            return dist

def exact_solution(X,t):
    radius = 0.15
    if nd==2:
        xc = 0.5
        yc = 0.75
        if ct.test_case==1:
            r = np.sqrt((X[0]-xc)**2 + (X[1]-yc)**2)
            return radius - r
        else: #Zalesak disk
            if type(X) != dict:
                return zalesak_disk_per_point(X,0)
            else:
                x = X[0]
                y = X[1]
                z=np.zeros(x.shape,'d')
                for i in range(len(x)):
                    j=0
                    for xq,yq in zip(x[i],y[i]):
                        XX = {0:xq,1:yq}
                        z[i,j] = zalesak_disk_per_point(XX,t)
                        j = j+1
                return z
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
class init_cond(object):
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return exact_solution(x,0)
initialConditions  = {0:init_cond()}

#######################
# BOUNDARY CONDITIONS #
#######################
def getDBC(x,flag):
    return lambda x,t: 0.
dirichletConditions = {0:getDBC}

def zeroadv(x,flag):
    None
advectiveFluxBoundaryConditions =  {0:zeroadv}

fluxBoundaryConditions = {0:'outFlow'}
diffusiveFluxBoundaryConditions = {0:{}}
