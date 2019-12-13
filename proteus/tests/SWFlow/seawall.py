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
import proteus.SWFlow.SWFlowProblem as SWFlowProblem
from proteus import WaveTools as wt


# This is the problem of a solitary wave overtopping
# a seawall. Experiments were conducted by
# Hsiao and Lin is 2010 (?) and set up of the problem
# can be found at http://mail178.taseptrev.com/src/test/seawall.c
# Note that this is a fake 1D problem, ie
# we are doing simulation in 2d but only consider x direction velocity

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts = Context.Options([
    ('sw_model', 1, "sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time", 12.0, "Final time for simulation"),
    ("dt_output", 0.1, "Time interval to output solution"),
    ("cfl", 0.25, "Desired CFL restriction"),
    ("refinement", 4, "Refinement level"),
    ("reflecting_BCs",False,"Use reflecting BCs")
])

###################
# DOMAIN AND MESH #
###################
L = (15.0, 1.0)  # this is length in x direction and y direction
refinement = opts.refinement
domain = RectangularDomain(L=L, x=[0, 0, 0])  # x is bottom left corner

# CREATE REFINEMENT #
nnx0 = 6
nnx = (nnx0 - 1) * (2**refinement) + 1
nny = old_div((nnx - 1), 10) + 1
he = old_div(L[0], float(nnx - 1))
triangleOptions="pAq30Dena%f"  % (0.5*he**2,)

###############################
#  CONSTANTS NEEDED FOR SETUP #
###############################
g = 9.81
# water depth and relative amplitude
h0 = 0.2
a = 0.35
# wavenumber, width and wave celerity
k_wavenumber = np.sqrt(3.0 * a / (4.0 * h0**3))
z = np.sqrt(3.0 * a * h0) / (2.0 * h0 * np.sqrt(h0 * (1.0 + a)))
c = np.sqrt(g * (1.0 + a) * h0)
# initial location of solitary wave
x0 = 5.9

###############################
#   Functions defined here    #
###############################


def solitary_wave(x, t):
    sechSqd = (1.00 / np.cosh(z * (x - x0 - c * t)))**2.00
    return a * h0 * sechSqd


def bathymetry_function(X):
    #  need this shift for experimental data
    x = X[0] + 3
    conds = [x < 10, (13.6 < x) & (x <= 13.9), (13.9 < x) & (x <= 13.948), \
                        (13.948 < x) & (x<= 14.045)]
    # a bit frustrating you can't do bath = bath - h0 so have to subtract h0 in
    # every function
    bath = [lambda x: 0 - h0, \
    lambda x: 3.6 / 20. + 0.076 / (13.9 - 13.6) * (x - 13.6) - h0, \
    lambda x: 3.6 / 20. + 0.076 - h0, \
    lambda x: 3.6 / 20. + 0.076 - (0.076 - 0.022) / (14.045 - 13.948) * (x - 13.948) - h0, \
    lambda x: 1 / 20. * (x - 10.) - h0]
    return np.piecewise(x, conds, bath)

##############################
##### INITIAL CONDITIONS #####
##############################
class water_height_at_t0(object):
    def uOfXT(self, X, t):
        eta = solitary_wave(X[0], 0)
        h = max(eta - bathymetry_function(X), 0.)
        return h


class x_mom_at_t0(object):
    def uOfXT(self, X, t):
        eta = solitary_wave(X[0], 0)
        h = max(eta - bathymetry_function(X), 0.)
        Umom = h * c * eta / (h0 + eta)
        return Umom


class y_mom_at_t0(object):
    def uOfXT(self, X, t):
        return 0.


"""
heta and hw are needed for the modified green naghdi equations.
Note that the BCs for the heta and hw should be same as h.
For more details see: 'Robust explicit relaxation techinque for solving
the Green-Naghdi equations' by Guermond, Popov, Tovar, Kees.
JCP 2019
"""

class heta_at_t0(object):
    def uOfXT(self, X, t):
        h = water_height_at_t0().uOfXT(X, t)
        return h**2


class hw_at_t0(object):
    def uOfXT(self, X, t):
        eta = solitary_wave(X[0], 0)
        h = max(eta - bathymetry_function(X), 0.)
        hprime = -2.0 * z * eta * np.tanh(z * (X[0] - x0 - c * t))
        hw = h * (-c * h0 * eta * hprime / (h0 + eta)**2)
        return hw

###############################
##### BOUNDARY CONDITIONS #####
###############################
X_coords = (0.0, 15.0)  # this is domain in x direction, used for BCs

def x_mom_DBC(X, flag):
    if X[0] == X_coords[0] or X[0] == X_coords[1]:
        return lambda x, t: 0.0

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
boundaryConditions = {'water_height': lambda x, flag: None,
                      'x_mom': x_mom_DBC,
                      'y_mom': lambda x, flag: lambda x, t: 0.0,
                      'h_times_eta': lambda x, flag: None,
                      'h_times_w': lambda x, flag: None}
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
mySWFlowProblem.physical_parameters['mannings'] = 0.012
