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
import proteus.SWFlows.SWFlowProblem as SWFlowProblem 

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ('sw_model',0,"sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time",1000.0,"Final time for simulation"),
    ("dt_output",10.0,"Time interval to output solution"),
    ("cfl",0.33,"Desired CFL restriction")
    ])

###################
# DOMAIN AND MESH #
###################
L=(10000.0,10000.0)
refinement = 3
domain = RectangularDomain(L=L,x=[0,0,0])

# CREATE REFINEMENT #
nnx0=6
nnx = nny = (nnx0-1)*(2**refinement)+1
he = old_div(L[0],float(nnx-1))
triangleOptions="pAq30Dena%f"  % (0.5*he**2,)

######################
##### BATHYMETRY #####
######################
h0=10
a=3000
B=5
k=0.002
g = SWFlowProblem.default_physical_parameters['gravity']
p = old_div(np.sqrt(8*g*h0),a)
s = old_div(np.sqrt(p**2 - k**2),2.)
mannings = k

def bathymetry_function(X):
    x = X[0]
    y = X[1] 
    r2 = (x-old_div(L[0],2.))**2+(y-old_div(L[1],2.))**2
    return h0*r2/a/a    

##############################
##### INITIAL CONDITIONS #####
##############################
def eta_function(X,t):
    x = X[0]
    y = X[1]

    coeff1 = old_div(-B,g)
    coeff2 = -B**2/2./g
    
    eta_part1 = coeff1*np.exp(-k*t/2.)*(k/2.*np.sin(s*t)+s*np.cos(s*t))*(x-old_div(L[0],2.));
    eta_part2 = coeff1*np.exp(-k*t/2.)*(k/2.*np.cos(s*t)-s*np.sin(s*t))*(y-old_div(L[1],2.));
    eta_part3 = coeff2*np.exp(-k*t);
    
    return h0 + eta_part1 + eta_part2 + eta_part3

class water_height_at_t0(object):
    def uOfXT(self,X,t):
        return max(eta_function(X,0)-bathymetry_function(X),0)
    
class x_mom(object):
    def uOfXT(self,X,t):
        return 0.

class y_mom(object):
    def uOfXT(self,X,t):
        h = water_height_at_t0().uOfXT(X,t)
        v = B
        return h*v    

# ********************************** #
# ***** Create mySWFlowProblem ***** #
# ********************************** #
outputStepping = SWFlowProblem.OutputStepping(opts.final_time,dt_output=opts.dt_output)
initialConditions = {'water_height': water_height_at_t0(),
                     'x_mom': x_mom(),
                     'y_mom': y_mom()}
boundaryConditions = {'water_height': lambda x,flag: None,
                      'x_mom': lambda x,flag: None,
                      'y_mom': lambda x,flag: None}
mySWFlowProblem = SWFlowProblem.SWFlowProblem(sw_model=0,
                                              cfl=0.33,
                                              outputStepping=outputStepping,
                                              structured=True,
                                              he=he,
                                              nnx=nnx,
                                              nny=nny,
                                              domain=domain,
                                              initialConditions=initialConditions,
                                              boundaryConditions=boundaryConditions,
                                              bathymetry=bathymetry_function)
mySWFlowProblem.physical_parameters['LINEAR_FRICTION']=1
mySWFlowProblem.physical_parameters['mannings']=mannings
