from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2DCV
from proteus.Domain import RectangularDomain
import numpy as np
from proteus import (Domain, Context,
                     MeshTools as mt)
from proteus.Profiling import logEvent
import proteus.SWFlows.SWFlowProblem as SWFlowProblem
from proteus import WaveTools as wt

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts = Context.Options([
    ('sw_model', 0, "sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time", 10.0, "Final time for simulation"),
    ("dt_output", 0.1, "Time interval to output solution"),
    ("cfl", 0.2, "Desired CFL restriction"),
    ("refinement", 4, "Refinement level")
])

###################
# DOMAIN AND MESH #
###################
L = (600.0, 1.0)  # this is length in x direction and y direction
refinement = opts.refinement
# x is bottom left corner of rectangle
domain = RectangularDomain(L=L, x=[-300, 0, 0])
X_coords = (-300.0, 300.0)  # this is domain in x direction, used for BCs

# CREATE REFINEMENT #
nnx0 = 6
nnx = (nnx0 - 1) * (2**refinement) + 1
nny = old_div((nnx - 1), 10) + 1

he = old_div(L[0], float(nnx - 1))
triangleOptions = "pAq30Dena%f" % (0.5 * he**2,)

###############################
#  CONSTANTS NEEDED FOR SETUP #
###############################
hL = 1.8  # hL for hLeft
hR = 1.0  # hR for hRight
g = 9.81


def bathymetry_function(X):
    return 0.*X[0]

##############################
##### INITIAL CONDITIONS #####
##############################


class water_height_at_t0(object):
    def uOfXT(self, X, t):
        x = X[0]
        if x < 0:
            h = hL
        else:
            h = hR
        return h


class x_mom_at_t0(object):
    def uOfXT(self, X, t):
        return 0.


class y_mom_at_t0(object):
    def uOfXT(self, X, t):
        return 0.* X[0]
# heta and hw are needed for the dispersive modified green naghdi equations
# source is 'ROBUST EXPLICIT RELAXATION TECHNIQUE FOR SOLVING
# THE GREEN NAGHDI EQUATIONS' by Guermond, Popov, Tovar
class heta_at_t0(object):
    def uOfXT(self, X, t):
        h = water_height_at_t0().uOfXT(X, t)
        return h**2
class hw_at_t0(object):
    def uOfXT(self, X, t):
        return 0.

###############################
##### BOUNDARY CONDITIONS #####
###############################
def water_height_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda x, t: water_height_at_t0().uOfXT(X, 0.0)
    elif X[0] == X_coords[1]:
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
boundaryConditions = {'water_height': water_height_DBC,
                      'x_mom': x_mom_DBC,
                      'y_mom': lambda x, flag: lambda x, t: 0.0,
                      'h_times_eta': heta_DBC,
                      'h_times_w': hw_DBC}
mySWFlowProblem = SWFlowProblem.SWFlowProblem(sw_model=0,
                                              cfl=0.2,
                                              outputStepping=outputStepping,
                                              structured=True,
                                              he=he,
                                              nnx=nnx,
                                              nny=nny,
                                              domain=domain,
                                              initialConditions=initialConditions,
                                              boundaryConditions=boundaryConditions,
                                              bathymetry=bathymetry_function)
mySWFlowProblem.physical_parameters['LINEAR_FRICTION'] = 0
mySWFlowProblem.physical_parameters['mannings'] = 0
