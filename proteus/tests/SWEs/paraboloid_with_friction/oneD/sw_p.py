from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2DCV
from proteus.Domain import RectangularDomain
import numpy as np

# ************************************** #
# ********** PROBLEM SPECIFIC ********** #
# ************************************** #
refinement=7

#domain and parameters
L=(8000.0,800.0)
g = 9.81
h0=10
a=3000
B=2
k=0.001
p = old_div(np.sqrt(8*g*h0),a)
s = old_div(np.sqrt(p**2 - k**2),2.)

# For time and outputting
T=1000.00
nDTout=10

AUTOMATED_TEST=True
if AUTOMATED_TEST:
    T=500.0
    nDTout=1
    refinement=3
    
######################
##### PARAMETERS #####
######################
# PHYSICAL PARAMETERS #
g = 9.81
LINEAR_FRICTION = 1
mannings = k

# NUMERICAL PARAMETERS #
cE = 1
LUMPED_MASS_MATRIX = 0
SSPOrder = 3
runCFL=0.25
useSuperlu = True
triangleFlag = 1
reflecting_BCs = False

##################
##### DOMAIN #####
##################
nd=2
domain = RectangularDomain(L=L,x=[0,0,0])
bt = domain.boundaryTags
bt['front'] = bt['bottom']
bt['back'] = bt['top']
domain.writePoly("tank2d")

################
##### MESH #####
################
nnx0=6
nnx = (nnx0-1)*(2**refinement)+1
nny = old_div((nnx-1),10)+1
nnz=1
he = old_div(L[0],float(nnx-1))
triangleOptions="pAq30Dena%f"  % (0.5*he**2,)
domain.MeshOptions.triangleOptions=triangleOptions

######################
##### BATHYMETRY #####
######################
def bathymetry(X):
    x=X[0]
    return h0*(x-old_div(L[0],2))**2/a/a

def eta_function(x,t):
    coeff1 = a**2*B**2/8./g/g/h0
    coeff2 = -B**2/4./g
    coeff3 = old_div(-1.,g)
    
    eta_part1 = coeff1*np.exp(-k*t)*(-s*k*np.sin(2*s*t)+(old_div(k**2,4.)-s**2)*np.cos(2*s*t))
    eta_part2 = coeff2*np.exp(-k*t)
    eta_part3 = coeff3*np.exp(-k*t/2.)*(B*s*np.cos(s*t)+k*B/2.*np.sin(s*t))*(x-old_div(L[0],2))
    
    return h0 + eta_part1 + eta_part2 + eta_part3

##############################
##### INITIAL CONDITIONS #####
##############################
class IC_h(object):
    def uOfXT(self,X,t):
        eta = eta_function(X[0],0)
        h = eta-bathymetry(X)
        hp = max(h,0.)
        return hp

class Zero(object):
    def uOfXT(self,x,t):
        return 0.0

IC_hu=Zero
IC_hv=Zero

##########################
##### EXACT SOLUTION #####
##########################
class water_height(object):
    def uOfXT(self,X,t):
        h = eta_function(X[0],t) - bathymetry(X)
        hp = max(h,0.)
        return hp
        
class exact_velx(object):
    def uOfXT(self,X,t):
        return B*np.exp(-k*t/2.)*np.sin(s*t)

class exact_momx(object):
    def uOfXT(self,X,t):
        ht = eta_function(X[0],t) - bathymetry(X)
        if (ht >= 0):
            h=ht
        else:
            h=0
        u = B*np.exp(-k*t/2.)*np.sin(s*t)
        return h*u

analyticalSolution = {0:water_height(),
                      1:exact_momx(),
                      2:Zero()}

###################################
##### FOR BOUNDARY CONDITIONS #####
###################################
def getDBC_h(x,flag):
    return None

def getDBC_hu(x,flag):
    return None
    
def getDBC_hv(x,flag):
    return lambda x,t: 0.0

# ************************************************ #
# ********** GENERIC PORTION OF _p FILE ********** #
# ************************************************ #

# ********************************** #
# ********** COEFFICIENTS ********** #
# ********************************** #
LevelModelType = SW2DCV.LevelModel
coefficients = SW2DCV.Coefficients(g=g,
                                   bathymetry={0:bathymetry},
                                   cE=cE,
                                   LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX,
                                   LINEAR_FRICTION=LINEAR_FRICTION,
                                   mannings=mannings)


# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions = {0: IC_h(),
                     1: IC_hu(),
                     2: IC_hv()}

# ***************************************** #
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
dirichletConditions = {0:getDBC_h,
                       1:getDBC_hu,
                       2:getDBC_hv}
fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}
advectiveFluxBoundaryConditions =  {0: lambda x,flag: None,
                                    1: lambda x,flag: None,
                                    2: lambda x,flag: None }
diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1: lambda x,flag: None},
                                   2:{2: lambda x,flag: None}}
