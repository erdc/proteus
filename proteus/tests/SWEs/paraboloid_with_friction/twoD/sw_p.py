from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2DCV
from proteus.Domain import RectangularDomain
import numpy as np
import math

# ************************************** #
# ********** PROBLEM SPECIFIC ********** #
# ************************************** #
refinement=6

#domain and parameters
L=(10000.0,10000.0)
g=9.81
h0=10.0
a=3000.
B=5.
k=0.002
p = old_div(np.sqrt(8*g*h0),a)
s = old_div(np.sqrt(p**2 - k**2),2.)
mannings=k

# For time and outputting
T=1000. 
nDTout=10

AUTOMATED_TEST=True
if AUTOMATED_TEST:
    T=500
    nDTout=1
    refinement=2
    
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
domain = RectangularDomain(L=L)
bt = domain.boundaryTags
bt['front'] = bt['bottom']
bt['back'] = bt['top']
domain.writePoly("tank2d")

################
##### MESH #####
################
nnx0=6
nnx = (nnx0-1)*(2**refinement)+1
nny = nnx
nnz = 1
he = old_div(L[0],float(nnx-1))
triangleOptions="pAq30Dena%f"  % (0.5*he**2,)
domain.MeshOptions.triangleOptions=triangleOptions

######################
##### BATHYMETRY #####
######################
def bathymetry(X):
    import numpy as np
    x = X[0]
    y = X[1] 
    r2 = (x-old_div(L[0],2.))**2+(y-old_div(L[1],2.))**2
    return h0*r2/a/a

# AUXILIARY FUNCTIONS #
def eta_function(X,t):
    x = X[0]
    y = X[1]

    coeff1 = old_div(-B,g)
    coeff2 = -B**2/2./g
    
    eta_part1 = coeff1*np.exp(-k*t/2.)*(k/2.*np.sin(s*t)+s*np.cos(s*t))*(x-old_div(L[0],2.));
    eta_part2 = coeff1*np.exp(-k*t/2.)*(k/2.*np.cos(s*t)-s*np.sin(s*t))*(y-old_div(L[1],2.));
    eta_part3 = coeff2*np.exp(-k*t);
    
    return h0 + eta_part1 + eta_part2 + eta_part3

def water_height_function(X,t):
    eta = eta_function(X,t)
    return max(eta-bathymetry(X),0)

##############################
##### INITIAL CONDITIONS #####
##############################
class IC_h():
    def uOfXT(self,X,t):
        return water_height_function(X,0)

class IC_hu(object):
    def uOfXT(self,X,t):
        return 0.

class IC_hv(object):
    def uOfXT(self,X,t):
        h = water_height_function(X,0)
        v = B
        return h*v

##########################
##### EXACT SOLUTION #####
##########################
class water_height(object):
    def uOfXT(self,X,t):
        return water_height_function(X,t)

class Zero(object):
    def uOfXT(self,X,t):
        return 0.

analyticalSolution = {0:water_height(),
                      1:Zero(), 
                      2:Zero()}

###################################
##### FOR BOUNDARY CONDITIONS #####
###################################
def getDBC_h(x,flag):
    return None

def getDBC_hu(x,flag):
    return None

def getDBC_hv(x,flag):
    return None

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
