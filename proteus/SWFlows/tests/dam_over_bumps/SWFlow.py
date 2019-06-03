from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2D
from proteus.mprans import SW2DCV
from proteus.mprans import GN_SW2DCV
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
    ("final_time",30.0,"Final time for simulation"),
    ("dt_output",0.1,"Time interval to output solution"),
    ("refinement",4,"Level of refinement"),
    ("cfl",0.33,"Desired CFL restriction"),
    ("reflecting_BCs",True,"Use reflecting BCs")
    ])

###################
# DOMAIN AND MESH #
###################
L=(75.0,30.0)
refinement = opts.refinement
domain = RectangularDomain(L=L)

# CREATE REFINEMENT #
nnx0=6
nnx = (nnx0-1)*(2**refinement)+1
nny = old_div((nnx-1),2)+1
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
    bump1 = 1-1./8*np.sqrt((x-30)**2+(y-6)**2)
    bump2 = 1-1./8*np.sqrt((x-30)**2+(y-24)**2)
    bump3 = 3-3./10*np.sqrt((x-47.5)**2+(y-15)**2)
    return np.maximum(np.maximum(np.maximum(0.,bump1),bump2),bump3)

##############################
##### INITIAL CONDITIONS #####
##############################
class water_height_at_t0(object):
    def uOfXT(self,X,t):
        x = X[0]
        if (x <= 16):
            eta=1.875
        else:
            eta=0.

        z = bathymetry_function(X)
        return max(eta - z,0.)

class Zero(object):
    def uOfXT(self,x,t):
        return 0.0

# heta and hw are needed for the modified Green-Naghdi equations,
# ie dispersive shallow water equations


class heta_at_t0(object):
    def uOfXT(self, X, t):
        h = water_height_at_t0().uOfXT(X, t)
        return h**2


class hw_at_t0(object):
    def uOfXT(self, X, t):
        eta = solitary_wave(X[0], 0)
        h = eta - bathymetry_function(X)
        hp = max(h, 0.)
        hprime = -2.0 * z * eta * np.tanh(z * (X[0] - x0 - c * t))
        hw = hp * (-c * h0 * eta * hprime / (h0 + eta)**2)
        return hw

# ********************************** #
# ***** Create mySWFlowProblem ***** #
# ********************************** #
outputStepping = SWFlowProblem.OutputStepping(opts.final_time,dt_output=opts.dt_output)
initialConditions = {'water_height': water_height_at_t0(),
                     'x_mom': Zero(),
                     'y_mom': Zero(),
                     'h_times_eta': heta_at_t0(),
                     'h_times_w': Zero()}
boundaryConditions = {'water_height': lambda x,flag: None,
                      'x_mom': lambda x,flag: None,
                      'y_mom': lambda x,flag: None,
                      'h_times_eta': lambda x,flag: None,
                      'h_times_w': lambda x,flag: None}
mySWFlowProblem = SWFlowProblem.SWFlowProblem(sw_model=opts.sw_model,
                                              cfl=0.33,
                                              outputStepping=outputStepping,
                                              structured=True,
                                              he=he,
                                              nnx=nnx,
                                              nny=nny,
                                              domain=domain,
                                              initialConditions=initialConditions,
                                              boundaryConditions=boundaryConditions,
                                              reflectingBCs=opts.reflecting_BCs,
                                              bathymetry=bathymetry_function)
mySWFlowProblem.physical_parameters['LINEAR_FRICTION']=0
mySWFlowProblem.physical_parameters['mannings']=0.02
