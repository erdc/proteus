from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2DCV
from proteus.mprans import GN_SW2DCV
from proteus.Domain import RectangularDomain
import numpy as np
from proteus import (Domain, Context,
                     MeshTools as mt)
from proteus.Profiling import logEvent
import proteus.SWFlows.SWFlowProblem as SWFlowProblem

"""
This is the problem of a solitary wave overtopping a canonical island.
This experiment was done at (insert ref here). There were several gauges
places in the tank measuring the water height
"""

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts = Context.Options([
    ('sw_model', 0, "sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time", 9.0, "Final time for simulation"),
    ("dt_output", 0.1, "Time interval to output solution"),
    ("cfl", 0.30, "Desired CFL restriction"),
    ("refinement", 4, "Refinement level"),
    ("reflecting_BCs",True,"Use reflecting BCs")
])

###################
# DOMAIN AND MESH #
###################
L = (30.0, 25.0)  # this is length in x direction and y direction
refinement = opts.refinement
domain = RectangularDomain(L=L)

# CREATE REFINEMENT #
nnx0 = 6
nnx = (nnx0 - 1) * (2**refinement) + 1
nny = old_div((nnx - 1), 2) + 1
he = old_div(L[0], float(nnx - 1))
triangleOptions="pAq30Dena%f"  % (0.5*he**2,)

###############################
#  CONSTANTS NEEDED FOR SETUP #
###############################
g = 9.81

# stuff for solitary wave
h0 = 0.32
alpha = 1 * h0
xs = 15.0 - 12.96
vel = np.sqrt(old_div(3 * alpha, (4 * h0**3)))

# stuff for cone bathymetry
htop = 0.625
rcone = 3.6
scone = 4.0
hcone = 0.9


###############################
#   Functions defined here    #
###############################
def solitary_wave(x, t):
    sechSqd = (1.00 / np.cosh(vel * (x - xs)))**2.00
    soliton = alpha * sechSqd
    return soliton


def bathymetry_function(X):
    x = X[0]
    y = X[1]
    radius = np.sqrt((x - 15.0)**2 + (y - 13.0)**2)

    # need to do this annoying thing for piecewise functions
    conds = [radius <= rcone, radius > rcone]
    bath = [lambda radius: np.minimum(htop, hcone-radius/scone), lambda radius: 0.0]
    return np.piecewise(radius, conds, bath)

##############################
##### INITIAL CONDITIONS #####
##############################
class water_height_at_t0(object):
    def uOfXT(self, X, t):
        eta = h0 + solitary_wave(X[0], 0)
        h = max(eta - bathymetry_function(X), 0.)
        return h


class x_mom_at_t0(object):
    def uOfXT(self, X, t):
        return solitary_wave(X[0], 0) * np.sqrt(g / h0)


class y_mom_at_t0(object):
    def uOfXT(self, X, t):
        return 0.


# heta and hw are needed for the modified green naghdi equations
# source is 'ROBUST EXPLICIT RELAXATION TECHNIQUE FOR SOLVING
# THE GREEN NAGHDI EQUATIONS' by Guermond, Popov, Tovar, Kees

class heta_at_t0(object):
    def uOfXT(self, X, t):
        h = water_height_at_t0().uOfXT(X, t)
        return h**2


class hw_at_t0(object):
    def uOfXT(self, X, t):
        hw = 2.0 * solitary_wave(X[0], 0) * vel * np.sinh(vel * (X[0]-xs))
        return hw

###############################
##### BOUNDARY CONDITIONS #####
###############################

# Actually don't need any of these
def water_height_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda x, t: water_height_at_t0().uOfXT(X, 0.0)

def x_mom_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda X, t: x_mom_at_t0().uOfXT(X, 0.0)


def y_mom_DBC(X, flag):
    return lambda x, t: 0.0


def heta_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda x, t: heta_at_t0().uOfXT(X, 0.0)


def hw_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda x, t: hw_at_t0().uOfXT(X, 0.0)


# ********************************** #
# ***** Create mySWFlowProblem ***** #
# ********************************** #
outputStepping = SWFlowProblem.OutputStepping(
    opts.final_time, dt_output=opts.dt_output)
initialConditions = {'water_height': water_height_at_t0(),
                     'x_mom': x_mom_at_t0(),
                     'y_mom': y_mom_at_t0(),
                     'h_times_eta': heta_at_t0(),
                     'h_times_w': hw_at_t0()}
boundaryConditions = {'water_height': lambda x,flag: None,
                      'x_mom': lambda x,flag: None,
                      'y_mom': lambda x,flag: None,
                      'h_times_eta': lambda x,flag: None,
                      'h_times_w': lambda x,flag: None}
mySWFlowProblem = SWFlowProblem.SWFlowProblem(sw_model=opts.sw_model,
                                              cfl=opts.cfl,
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
mySWFlowProblem.physical_parameters['LINEAR_FRICTION'] = 0
mySWFlowProblem.physical_parameters['mannings'] = 0.0
