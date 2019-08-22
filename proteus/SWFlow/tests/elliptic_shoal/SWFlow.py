from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2DCV
from proteus.mprans import DSW2DCV
from proteus.Domain import RectangularDomain
import numpy as np
import math
from proteus import (Domain, Context,
                     MeshTools as mt,
                     WaveTools as wt)
from proteus.Profiling import logEvent
import proteus.SWFlow.SWFlowProblem as SWFlowProblem

"""
This is the problem of monochromatic waves propagating
over an elliptic shoal. Details can be found in section 4.5 of
A discontinuous Galerkin method for a new class of Green Naghdi
equations on simplicial unstructured meshes by Duran and Marche.
We use proteus's wave tools feature for defining
the monochromatic waves at the inlet.
"""

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts = Context.Options([
    ('sw_model', 0, "sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time", 10.0, "Final time for simulation"),
    ("dt_output", 0.1, "Time interval to output solution"),
    ("cfl", 0.33, "Desired CFL restriction"),
    ("refinement", 5, "Refinement level")
])

###################
# DOMAIN AND MESH #
###################
L = (15.0, 1.0)
domain = RectangularDomain(L=L, x=[0, 0, 0])
X_coords = (0.0, 15.0)  # this is domain in x direction, used in BCs

# CREATE REFINEMENT #
refinement = opts.refinement
nnx0 = 6
nnx = (nnx0 - 1) * (2**refinement) + 1
nny = old_div((nnx - 1), 10) + 1
he = old_div(L[0], float(nnx - 1))
triangleOptions = "pAq30Dena%f" % (0.5 * he**2,)

###############################
#  Wave tools set up + bath   #
###############################
h0 = 0.38  # depth / still water level
amp = 0.0232  # amplitude of waves / wave height
period = 1  # period in seconds
wave_dir = (1., 0., 0.)  # wave direction
g_vec = [0, -9.81, 0]  # g vector used for wave tools

# define wave here
wave = wt.MonochromaticWaves(
    period, amp, h0, h0, np.array(g_vec), wave_dir, "Fenton")


def bathymetry_function(X):
    slope = 1.0 / 20.
    x = X[0]
    bath = 0.5 - 0.2 + 0.0 * x
    # silly hack because X switches from list to array of
    # length 3 (x,y,z) when called in initial conditions
    if (isinstance(X, list)):
        for counter, value in enumerate(X[0]):
            if value < 10.0:
                bath[counter] = np.maximum(slope * value, 0.2) - 0.2
    else:
        if (x <= 10.):
            bath = np.maximum(slope * x, 0.2) - 0.2
    return bath


##############################
##### INITIAL CONDITIONS #####
##############################
class water_height_at_t0(object):
    def uOfXT(self, X, t):
        h = np.maximum(h0 - bathymetry_function(X), 0.0)
        return h


class x_mom_at_t0(object):
    def uOfXT(self, X, t):
        return 0.0


class y_mom_at_t0(object):
    def uOfXT(self, X, t):
        return 0.0


# heta and hw are needed for the dispersive modified green naghdi equations
# source is 'ROBUST EXPLICIT RELAXATION TECHNIQUE FOR SOLVING
# THE GREEN NAGHDI EQUATIONS' by Guermond, Popov, Tovar


class heta_at_t0(object):
    def uOfXT(self, X, t):
        return 0.0


class hw_at_t0(object):
    def uOfXT(self, X, t):
        return 0.0


###############################
##### BOUNDARY CONDITIONS #####
###############################


def water_height_DBC(X, flag):
    return None


def x_mom_DBC(X, flag):
    if X[0] == X_coords[0]:
        h = np.maximum(h0 - bathymetry_function(X), 0.0)
        return lambda X, t: h * wave.u(X, t)[0]
    elif X[0] == X_coords[1]:
        return lambda X, t: 0.0


def y_mom_DBC(X, flag):
    return lambda x, t: 0.0


def heta_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda x, t: heta_at_t0().uOfXT(X, 0.0)  # FIX FOR MGN, EJT


def hw_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda x, t: hw_at_t0().uOfXT(X, 0.0)  # FIX FOR MGN, EJT


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
                                              bathymetry=bathymetry_function)
mySWFlowProblem.physical_parameters['LINEAR_FRICTION'] = 0
mySWFlowProblem.physical_parameters['mannings'] = 0.0
