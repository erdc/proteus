from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2D
from proteus.mprans import SW2DCV
from proteus.Domain import RectangularDomain
import numpy as np
import math

nd=2

T=10.#2*math.pi/omega
nDTout=100

L=(10.0,10.0)
g = 9.81
# PARAMETERS #
h0=1.0
a=3
B=2
k=0.5
p = np.sqrt(8*g*h0)/a
s = np.sqrt(p**2 - k**2)/2.
mannings=k

domain = RectangularDomain(L=L)

#This is relevant just when use_second_order_NonFlatB_with_EV_stabilization=True
cE=1
LUMPED_MASS_MATRIX=1

bt = domain.boundaryTags
bt['front'] = bt['bottom']
bt['back'] = bt['top']
domain.writePoly("tank2d")

######################
##### BATHYMETRY #####
######################
def bathymetry_function(X):
    import numpy as np
    x = X[0]
    y = X[1] 
    r2 = (x-L[0]/2.)**2+(y-L[0]/2.)**2
    return h0*r2/a/a

def eta_function(X,t):
    x = X[0]
    y = X[1]

    coeff1 = -B/g
    coeff2 = -B**2/2./g
    
    eta_part1 = coeff1*np.exp(-k*t/2.)*(k/2.*np.sin(s*t)+s*np.cos(s*t))*(x-L[0]/2.);
    eta_part2 = coeff1*np.exp(-k*t/2.)*(k/2.*np.cos(s*t)-s*np.sin(s*t))*(y-L[0]/2.);
    eta_part3 = coeff2*np.exp(-k*t);
    
    return h0 + eta_part1 + eta_part2 + eta_part3

def water_height_function(X,t):
    eta = eta_function(X,t)
    return max(eta-bathymetry_function(X),0)

##############################
##### INITIAL CONDITIONS #####
##############################
class water_height_at_t0:
    def uOfXT(self,X,t):
        return water_height_function(X,0)

class momX:
    def uOfXT(self,X,t):
        return 0.

class momY:
    def uOfXT(self,X,t):
        h = water_height_function(X,t)
        v = -B
        return h*v

class Zero:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:water_height_at_t0(),
                     1:momX(),
                     2:momY()}

##########################
##### EXACT SOLUTION #####
##########################
class water_height:
    def uOfXT(self,X,t):
        return water_height_function(X,t)

analyticalSolution = {0:water_height(),
                      1:Zero(), 
                      2:Zero()}

###################################
##### FOR BOUNDARY CONDITIONS #####
###################################
def getDBC_h(x,flag):
    return None

#note, these are the same for hu and hv so we can cheat and use  this p-file for SW2DCV and SW2D
def getDBC_u(x,flag):
    if (x[0] in [0.0,L[0]]) or flag in [bt['left'],bt['right']]:
        return lambda x,t: 0.0

def getDBC_v(x,flag):
    if x[1] in [0.0,L[1]] or flag in [bt['front'],bt['back']]:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_h,
                       1:getDBC_u,
                       2:getDBC_v}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

def getAFBC_h(x,flag):
    return lambda x,t: 0.0

def getAFBC_u(x,flag):
    if flag == 0:
        return lambda x,t: 0.0
    else:
        return None
def getAFBC_v(x,flag):
    if flag == 0:
        return lambda x,t: 0.0
    else:
        return None

advectiveFluxBoundaryConditions =  {0:getAFBC_h,
                                    1:getAFBC_u,
                                    2:getAFBC_v}

def getDFBC_u(x,flag):
    if flag == 0:
        return lambda x,t: 0.0
    else:
        return None

def getDFBC_v(x,flag):
    if flag == 0:
        return lambda x,t: 0.0
    else:
        return None

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v}}

#########################################
##### CREATE MODEL AND COEFFICIENTS #####
#########################################
bathymetry={0:bathymetry_function}
LevelModelType = SW2DCV.LevelModel
coefficients = SW2DCV.Coefficients(g=g,bathymetry=bathymetry,cE=cE,LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX,mannings=mannings)

