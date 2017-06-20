from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2D
from proteus.mprans import SW2DCV
from proteus.Domain import RectangularDomain
import numpy as np

nd=2

T=10
nDTout=100

L=(10.0,1.0)
g = 9.81
domain = RectangularDomain(L=L,x=[0,0,0])

h0=1
a=3
B=1
k=0.5
p = np.sqrt(8*g*h0)/a
s = np.sqrt(p**2 - k**2)/2.
mannings=k

#This is relevant just when use_second_order_NonFlatB_with_EV_stabilization=True
cE=1000000
LUMPED_MASS_MATRIX=0

bt = domain.boundaryTags
bt['front'] = bt['bottom']
bt['back'] = bt['top']
domain.writePoly("tank2d")

######################
##### BATHYMETRY #####
######################
def bathymetry_function(X):
    x=X[0]
    return h0*(x-L[0]/2)**2/a/a

def eta_function(x,t):
    coeff1 = a**2*B**2/8./g/g/h0
    coeff2 = -B**2/4./g
    coeff3 = -1./g
    
    eta_part1 = coeff1*np.exp(-k*t)*(-s*k*np.sin(2*s*t)+(k**2/4.-s**2)*np.cos(2*s*t))
    eta_part2 = coeff2*np.exp(-k*t)
    eta_part3 = coeff3*np.exp(-k*t/2.)*(B*s*np.cos(s*t)+k*B/2.*np.sin(s*t))*(x-L[0]/2)
    
    return h0 + eta_part1 + eta_part2 + eta_part3

##############################
##### INITIAL CONDITIONS #####
##############################
class water_height_at_t0:
    def uOfXT(self,X,t):
        eta = eta_function(X[0],0)
        h = eta-bathymetry_function(X)
        hp = max(h,0.)
        return hp

class Zero:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:water_height_at_t0(),
                     1:Zero(),
                     2:Zero()}

##########################
##### EXACT SOLUTION #####
##########################
class water_height:
    def uOfXT(self,X,t):
        h = eta_function(X[0],t) - bathymetry_function(X)
        hp = max(h,0.)
        return hp
        
class exact_velx:
    def uOfXT(self,X,t):
        return B*np.exp(-k*t/2.)*np.sin(s*t)

class exact_momx:
    def uOfXT(self,X,t):
        ht = eta_function(X[0],t) - bathymetry_function(X)
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
#    if (x[0] in [0.0]) or flag in [bt['left']]:
#        return lambda x,t: eta_function(x[0],t) - bathymetry_function(x)
#    else:
#        return None

def getDBC_u(x,flag):
    return None
#    if (x[0] in [0.0]) or flag in [bt['left']]:
#        return lambda x,t: (eta_function(x[0],t) - bathymetry_function(x)) * (B*np.exp(-k*t/2.)*np.sin(s*t))
#    else:
#        return None
    
def getDBC_v(x,flag):
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

