from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2D
from proteus.mprans import SW2DCV
from proteus.Domain import RectangularDomain
import numpy as np

nd=2

T=1000.00
nDTout=10

L=(8000.0,800.0)
g = 9.81
domain = RectangularDomain(L=L,x=[0,0,0])

h0=10
a=3000
B=2
k=0.001
p = old_div(np.sqrt(8*g*h0),a)
s = old_div(np.sqrt(p**2 - k**2),2.)
mannings=k

cE=1
LUMPED_MASS_MATRIX=0
LINEAR_FRICTION=1

bt = domain.boundaryTags
bt['front'] = bt['bottom']
bt['back'] = bt['top']
domain.writePoly("tank2d")

######################
##### BATHYMETRY #####
######################
def bathymetry_function(X):
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
class water_height_at_t0(object):
    def uOfXT(self,X,t):
        eta = eta_function(X[0],0)
        h = eta-bathymetry_function(X)
        hp = max(h,0.)
        return hp

class Zero(object):
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:water_height_at_t0(),
                     1:Zero(),
                     2:Zero()}

##########################
##### EXACT SOLUTION #####
##########################
class water_height(object):
    def uOfXT(self,X,t):
        h = eta_function(X[0],t) - bathymetry_function(X)
        hp = max(h,0.)
        return hp
        
class exact_velx(object):
    def uOfXT(self,X,t):
        return B*np.exp(-k*t/2.)*np.sin(s*t)

class exact_momx(object):
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
coefficients = SW2DCV.Coefficients(g=g,bathymetry=bathymetry,cE=cE,LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX,LINEAR_FRICTION=LINEAR_FRICTION,mannings=mannings)

