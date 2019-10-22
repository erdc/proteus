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
from proteus import WaveTools as wt

"""
This is the problem of a solitary wave run-up on a sloping beach.
Experiments were reported in Synolakas 1987, 'The runup of solitary waves'.
We use the probelm set up seen in http://basilisk.fr/src/test/beach.c
"""

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #

# Here we define Tstar which corresponds to experimental data time
# Tstar = {10,15,20,...,65} for the final time

T = 65.0
g = 9.81  # gravity
h0 = 1.0  # water depth
Tstar = T * np.sqrt(h0 / g)

opts = Context.Options([
    ('sw_model', 1, "sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time", Tstar, "Final time for simulation"),
    ("dt_output", 0.1, "Time interval to output solution"),
    ("cfl", 0.2, "Desired CFL restriction"),
    ("refinement", 4, "Refinement level"),
    ("reflecting_BCs",False,"Use reflecting BCs")
])

###################
# DOMAIN AND MESH #
###################
L = (50.0, 1.0)
refinement = opts.refinement
domain = RectangularDomain(L=L, x=[-35.0, 0, 0])
X_coords = (-35.0, 15.0)  # this is domain in x direction, used in BCs

# CREATE REFINEMENT #
nnx0 = 6
nnx = (nnx0 - 1) * (2**refinement) + 1
nny = old_div((nnx - 1), 10) + 1
he = old_div(L[0], float(nnx - 1))
triangleOptions = "pAq30Dena%f" % (0.5 * he**2,)

###############################
#  CONSTANTS NEEDED FOR SETUP #
###############################

a = 0.28  # relative amplitude
slope = 1.0 / 19.85  # beach slope
k_wavenumber = np.sqrt(3.0 * a / (4.0 * h0**3))  # wavenumber
z = np.sqrt(3.0 * a * h0) / (2.0 * h0 * np.sqrt(h0 * (1.0 + a)))  # width of solitary wave
L_wave = 2.0 / k_wavenumber * np.arccosh(np.sqrt(20.0))  # wavelength of solitary wave
c = np.sqrt(g * (1.0 + a) * h0)  # wave celerity (or speed)
x0 = - h0 / slope - L_wave / 2.0  # location of the toe of the beach

###############################
#   Functions defined here    #
###############################


def solitary_wave(x, t):
    sechSqd = (1.00 / np.cosh(z * (x - x0 - c * t)))**2.00
    return a * h0 * sechSqd


def bathymetry_function(X):
    x = X[0]
    return np.maximum(slope * x, -h0)


class Zero(object):
    def uOfXT(self, x, t):
        return 0.0


######################
# INITIAL CONDITIONS #
######################


class water_height_at_t0(object):
    def uOfXT(self, X, t):
        hTilde = solitary_wave(X[0], 0)
        h = max(hTilde - bathymetry_function(X), 0.0)
        return h


class x_mom_at_t0(object):
    def uOfXT(self, X, t):
        hTilde = solitary_wave(X[0], 0)
        h = max(hTilde - bathymetry_function(X), 0.0)
        return h * c * hTilde / (h0 + hTilde)


class y_mom_at_t0(object):
    def uOfXT(self, X, t):
        h = water_height_at_t0().uOfXT(X, t)
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
        hTilde = solitary_wave(X[0], 0)
        h = max(hTilde - bathymetry_function(X), 0.0)
        hTildePrime = -2.0 * z * hTilde * np.tanh(z * (X[0] - x0 - c * t))
        hw = -h**2 * (c * h0 * hTildePrime / (h0 + hTilde)**2)
        return hw

###############################
##### BOUNDARY CONDITIONS #####
###############################


def water_height_DBC(X, flag):
    return None

def x_mom_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda X, t: 0.0
    if X[0] == X_coords[1]:
        return lambda X, t: 0.0


def y_mom_DBC(X, flag):
    return lambda x, t: 0.0


def heta_DBC(X, flag):
    return None


def hw_DBC(X, flag):
    return None


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
#analyticalSolution
analyticalSolution={'h_exact': Zero(),
                  'hu_exact': Zero(),
                  'hv_exact': Zero(),
                  'heta_exact':Zero(),
                  'hw_exact':Zero()}
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
                                              bathymetry=bathymetry_function,
                                              analyticalSolution=analyticalSolution)
mySWFlowProblem.physical_parameters['LINEAR_FRICTION'] = 0
mySWFlowProblem.physical_parameters['mannings'] = 0.016
# mySWFlowProblem.swe_parameters['LUMPED_MASS_MATRIX'] = 1
