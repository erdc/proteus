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
numerical_parameters = {'LUMPED_MASS_MATRIX': 1,
                        'cfl': 0.33,
                        'SSPOrder': 3,
                        'cE': 10}
physical_parameters = {'g': 9.81,
                       'LINEAR_FRICTION': 0,
                       'mannings': 0.0}
###################
# OUTPUTTING TIME #
###################
T=4.00
nDTout=40

###################
# DOMAIN AND MESH #
###################
L=(10.0,1.0)
refinement = 4
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
    return 0.*X[0]

h1=0.10
h2=0.11
x0 = 2.0
D = np.sqrt(g * h2)

def solitary(X,t):
    xi = X - D*t-x0
    z1 = 3.0*(h2-h1)
    z2 = h2 * h1**2
    z = np.sqrt(z1 / z2)
    soliton =  h1 + (h2 - h1) * 1.0/(np.cosh(xi/2.0 * z)**2)
    return soliton
                  
##############################
##### INITIAL CONDITIONS #####
##############################
class water_height_at_t0(object):
    def uOfXT(self,X,t):
        eta = solitary(X[0],0)
        h = eta-bathymetry_function(X)
        hp = max(h,0.)
        return hp

class momentum_in_x_at_t0(object):
    def uOfXT(self,X,t):
        return D*(solitary(X[0],0) - h1)

class Zero(object):
    def uOfXT(self,X,t):
        return 0.0

initialConditions = {'water_height': water_height_at_t0(),
                     'x-mom': momentum_in_x_at_t0(),
                     'y-mom': Zero()}

def water_DBC(x,flag):
    if x[0]==0 or x[0]==L[0]:
        return lambda x,t: h1

def xmom(x,flag):
    if x[0]==0 or x[0]==L[0]:
        return lambda x,t: 0.0
    
###################################
##### FOR BOUNDARY CONDITIONS #####
###################################
boundaryConditions = {'water_height': water_DBC,
                      'x-mom': xmom,
                      'y-mom': lambda x,flag: None}
