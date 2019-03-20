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


"""
Note that this set up for a solitary wave is from  equation 31
in the the paper
A rapid numerical method for solving Serre Green Naghdi
equations describing long free surface gravity waves by
Favrie and Gavrilyuk. It assumes the reference height is h1.
"""

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #

opts = Context.Options([
    ('sw_model', 0, "sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time", 4.0, "Final time for simulation"),
    ("dt_output", 0.1, "Time interval to output solution"),
    ("cfl", 0.33, "Desired CFL restriction"),
    ("refinement", 4, "Refinement level")
])

###################
# DOMAIN AND MESH #
###################
L = (12.0, 1.0)
X_coords = (0.0, 12.0)  # this is domain, used in BCs
domain = RectangularDomain(L=L, x=[0, 0, 0])

# CREATE REFINEMENT #
refinement = opts.refinement
nnx0 = 6
nnx = (nnx0 - 1) * (2**refinement) + 1
nny = old_div((nnx - 1), 10) + 1
he = old_div(L[0], float(nnx - 1))
triangleOptions = "pAq30Dena%f" % (0.5 * he**2,)

###################################
# SOLITARY WAVE FUCTIONS AND BATH #
###################################
h1 = 0.1
h2 = 0.11
x0 = 2.0  # initial location of solitary wave
g = 9.81
D = np.sqrt(g * h2)
z = np.sqrt(old_div(3.0 * (h2 - h1), h2 * h1**2))


def soliton(x, t):
    phase = x - D * t - x0
    a1 = old_div(z * phase, 2.0)
    return h1 + (h2 - h1) * old_div(1.0, np.cosh(a1)**2)


def u(x, t):
    D = np.sqrt(g * h2)
    return D * (1.0 - old_div(h1, soliton(x, t)))


def bathymetry_function(X):
    x = X[0]
    return x * 0.0  # just 0


###################################
#    FOR ANALYTICAL SOLUTIONS     #
###################################
class Zero(object):
    def uOfXT(self, x, t):
        return 0.0


class water_height_at_tfinal(object):
    def uOfXT(self, X, t):
        return soliton(X[0], opts.final_time)


##############################
#    INITIAL CONDITIONS      #
##############################


class water_height_at_t0(object):
    def uOfXT(self, X, t):
        return soliton(X[0], 0.0)


class x_mom_at_t0(object):
    def uOfXT(self, X, t):
        return soliton(X[0], 0.0) * u(X[0], 0.0)


class y_mom_at_t0(object):
    def uOfXT(self, X, t):
        h = water_height_at_t0().uOfXT(X, t)
        return 0. * h


"""
heta and hw are needed for the dispersive modified green naghdi equations
source is 'ROBUST EXPLICIT RELAXATION TECHNIQUE FOR SOLVING
THE GREEN NAGHDI EQUATIONS' by Guermond, Kees, Popov, Tovar
"""


class heta_at_t0(object):
    def uOfXT(self, X, t):
        h = water_height_at_t0().uOfXT(X, 0.0)
        return h**2


class hw_at_t0(object):
    def uOfXT(self, X, t):
        # since there is no bathymatry hw = -h^2 * div(vel) = -D * h1 * h'
        x = X[0]
        D = np.sqrt(g * h2)
        z = np.sqrt(old_div(3.0 * (h2 - h1), h2 * h1**2))
        phase = x - D * t - x0
        a1 = old_div(z * phase, 2.0)
        sechSqd = old_div(1.0, np.cosh(a1)**2)
        hprime = -2.0 * (h2 - h1) * z * sechSqd * np.tanh(phase)
        hw = -D * h1 * hprime
        return hw


###############################
#     BOUNDARY CONDITIONS     #
###############################


def water_height_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda x, t: water_height_at_t0().uOfXT(X, 0.0)
    elif X[0] == X_coords[1]:
        return lambda x, t: water_height_at_t0().uOfXT(X, 0.0)


def x_mom_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda x, t: x_mom_at_t0().uOfXT(X, 0.0)
    elif X[0] == X_coords[1]:
        return lambda x, t: x_mom_at_t0().uOfXT(X, 0.0)


def y_mom_DBC(X, flag):
    return lambda x, t: 0.0


def heta_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda x, t: heta_at_t0().uOfXT(X, 0.0)
    elif X[0] == X_coords[1]:
        return lambda x, t: heta_at_t0().uOfXT(X, 0.0)


def hw_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda x, t: hw_at_t0().uOfXT(X, 0.0)
    elif X[0] == X_coords[1]:
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
if (opts.sw_model == 1):
    analyticalSolution = {0: water_height_at_tfinal(),
                          1: Zero(),
                          2: Zero(),
                          3: Zero(),
                          4: Zero()}
if (opts.sw_model == 0):
    analyticalSolution = {0: water_height_at_tfinal(),
                          1: Zero(),
                          2: Zero()}
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
                                              bathymetry=bathymetry_function,
                                              analyticalSolution=analyticalSolution)
mySWFlowProblem.physical_parameters['LINEAR_FRICTION'] = 0
mySWFlowProblem.physical_parameters['mannings'] = 0.0
