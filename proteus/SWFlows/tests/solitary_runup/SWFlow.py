from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2DCV
from proteus.mprans import DSW2DCV
from proteus.Domain import RectangularDomain
import numpy as np
from proteus import (Domain, Context,
                     MeshTools as mt)
from proteus.Profiling import logEvent
import proteus.SWFlows.SWFlowProblem as SWFlowProblem
from proteus import WaveTools as wt


# This is the problem of a solitary wave going up
# a plane beach. Experiments were conducted by
# Synolakas is 1980 (?) and set up of the problem
# can be found at http://mail178.taseptrev.com/src/test/beach.c
# Note that this is a fake 1D problem, ie
# we are doing simulation in 2d but only consider x direction velocity

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
T = 10.0
g = 9.81
h0 = 1.0
# Tstar corresponds to experimental data time.
Tstar = T * np.sqrt(h0 / g)

opts = Context.Options([
    ('sw_model', 0, "sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time", Tstar, "Final time for simulation"),
    ("dt_output", 0.1, "Time interval to output solution"),
    ("cfl", 0.1, "Desired CFL restriction"),
    ("refinement", 4, "Refinement level")
])
print opts.sw_model

###################
# DOMAIN AND MESH #
###################
L = (50.0, 1.0)
refinement = opts.refinement
domain = RectangularDomain(L=L, x=[-35.0, 0, 0])
X_coords = (-35.0, 15.0)  # this is domain, used in BCs

# CREATE REFINEMENT #
nnx0 = 6
nnx = (nnx0 - 1) * (2**refinement) + 1
nny = old_div((nnx - 1), 10) + 1

he = old_div(L[0], float(nnx - 1))
triangleOptions = "pAq30Dena%f" % (0.5 * he**2,)

###############################
#  CONSTANTS NEEDED FOR SETUP #
###############################

a = 0.30  # amplitude
slope = 1.0 / 19.850
k_wavenumber = np.sqrt(3.0 * a / (4.0 * h0**3))  # wavenumber
z = np.sqrt(3.0 * a * h0) / (2.0 * h0 * np.sqrt(h0 * (1.0 + a)))
# wavelength of solitary wave
L_wave = 2.0 / k_wavenumber * np.arccosh(np.sqrt(1.0 / 0.050))
c = np.sqrt(g * (1.0 + a) * h0)
x0 = - h0 / slope - L_wave / 2.0  # location of the toe of the beach

###############################
#   Functions defined here    #
###############################


def solitary_wave(x, t):
    sechSqd = (1.00 / np.cosh(z * (x - x0 - c * t)))**2.00
    soliton = a * h0 * sechSqd
    return soliton


def bathymetry_function(X):
    x = X[0]
    return numpy.maximum(slope * x, -h0)

##############################
##### INITIAL CONDITIONS #####
##############################


class water_height_at_t0(object):
    def uOfXT(self, X, t):
        eta = solitary_wave(X[0], 0)
        h = eta - bathymetry_function(X)
        hp = max(h, 0.)
        return hp


class x_mom_at_t0(object):
    def uOfXT(self, X, t):
        eta = solitary_wave(X[0], 0)
        h = eta - bathymetry_function(X)
        hp = max(h, 0.)
        Umom = hp * c * eta / (h0 + eta)
        return Umom


class y_mom_at_t0(object):
    def uOfXT(self, X, t):
        h = water_height_at_t0().uOfXT(X, t)
        return 0.
# heta and hw are needed for the dispersive modified green naghdi equations
# source is 'ROBUST EXPLICIT RELAXATION TECHNIQUE FOR SOLVING
# THE GREEN NAGHDI EQUATIONS' by Guermond, Popov, Tovar


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

###############################
##### BOUNDARY CONDITIONS #####
###############################


def water_height_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda x, t: water_height_at_t0().uOfXT(X, 0.0)
    # elif X[0]==X_coords[1]:
    #     return lambda x,t: water_height_at_t0().uOfXT(X,0.0)


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
mySWFlowProblem = SWFlowProblem.SWFlowProblem(sw_model=opts.sw_model,
                                              cfl=0.1,
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
mySWFlowProblem.physical_parameters['mannings'] = 0.016
