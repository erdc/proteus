from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2D
from proteus.mprans import SW2DCV
from proteus.Domain import RectangularDomain
import numpy as np
from proteus import (Domain, Context,
                     MeshTools as mt)
from proteus.Profiling import logEvent

#####################
# DEFINE PARAMETERS #
#####################
numerical_parameters = {'LUMPED_MASS_MATRIX': 0,
                        'cfl': 0.33,
                        'SSPOrder': 3,
                        'cE': 1}
physical_parameters = {'g': 9.81,
                       'LINEAR_FRICTION': 0,
                       'mannings': 0.0}
###################
# OUTPUTTING TIME #
###################
T=1000.00
nDTout=10

###################
# DOMAIN AND MESH #
###################
L=(8000.0,800.0)
refinement = 3
domain = RectangularDomain(L=L,x=[0,0,0])

# CREATE REFINEMENT #
nnx0=6
nnz=1
nnx = (nnx0-1)*(2**refinement)+1
nny = old_div((nnx-1),10)+1

he = old_div(L[0],float(nnx-1))
triangleOptions="pAq30Dena%f"  % (0.5*he**2,)

######################
##### BATHYMETRY #####
######################
h0=10
a=3000
B=2
k=0.001
g=physical_parameters['g']
p = old_div(np.sqrt(8*g*h0),a)
s = old_div(np.sqrt(p**2 - k**2),2.)

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
    def uOfXT(self,X,t):
        return 0.0

initialConditions = {'water_height': water_height_at_t0(),
                     'x-mom': Zero(),
                     'y-mom': Zero()}

###################################
##### FOR BOUNDARY CONDITIONS #####
###################################
boundaryConditions = {'water_height': lambda x,flag: None,
                      'x-mom': lambda x,flag: None,
                      'y-mom': lambda x,flag: lambda x,t: 0.0}
