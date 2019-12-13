from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2DCV
from proteus.mprans import GN_SW2DCV
from proteus.Domain import RectangularDomain
import numpy as np
from proteus import (Domain, Context, MeshTools as mt)
from proteus.Profiling import logEvent
import proteus.SWFlow.SWFlowProblem as SWFlowProblem


"""
This test uses the definition of a solitary wave from eqn 31
in the the paper 'A rapid numerical method for solving Serre Green Naghdi
equations describing long free surface gravity waves' by
Favrie and Gavrilyuk.
"""

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #

opts = Context.Options([
    ('sw_model', 1, "sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time", 4.0, "Final time for simulation"),
    ("dt_output", 0.1, "Time interval to output solution"),
    ("cfl", 0.25, "Desired CFL restriction"),
    ("refinement", 4, "Refinement level")
])

###################
# DOMAIN AND MESH #
###################
L = (10.0, 1.0)
X_coords = (0.0, 10.0)  # this is domain, used in BCs
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
g = 9.81
h1 = .10
h2 = .11
x0 = 2.0 # initial location of solitary wave
# solitary wave celerity and width
c = np.sqrt(g * h2)
r = np.sqrt(old_div(3.0 * (h2 - h1), 4 * h2 * h1**2))

def solitary_wave(x, t):
    phase = x - c * t - x0
    return h1 + (h2 - h1) * old_div(1.0, np.cosh(r * phase)**2)


def u(x, t):
    h = solitary_wave(x,t)
    return c * (1.0 - old_div(h1, h))


def bathymetry_function(X):
    x = X[0]
    # then return vector of zeros
    return x * 0.0


###################################
#    FOR ANALYTICAL SOLUTIONS     #
###################################
class Zero(object):
    def uOfXT(self, x, t):
        return 0.0


class water_height_at_tfinal(object):
    def uOfXT(self, X, t):
        return solitary_wave(X[0], opts.final_time)


##############################
#    INITIAL CONDITIONS      #
##############################


class water_height_at_t0(object):
    def uOfXT(self, X, t):
        return solitary_wave(X[0], 0.0)


class x_mom_at_t0(object):
    def uOfXT(self, X, t):
        h = solitary_wave(X[0], 0.0)
        return h * u(X[0], 0.0)


class y_mom_at_t0(object):
    def uOfXT(self, X, t):
        h = water_height_at_t0().uOfXT(X, t)
        return 0. * h


"""
heta and hw are needed for the modified green naghdi equations.
Note that the BCs for the heta and hw should be same as h.
For more details see: 'Robust explicit relaxation techinque for solving
the Green-Naghdi equations' by Guermond, Popov, Tovar, Kees.
JCP 2019
"""


class heta_at_t0(object):
    def uOfXT(self, X, t):
        h = water_height_at_t0().uOfXT(X, 0.0)
        return h**2


class hw_at_t0(object):
    def uOfXT(self, X, t):
        # since there is no bathymetry, waterHeight = htilde
        # hw = -waterHeight^2 * div(vel)
        #    = -h^2 * (c * h1 * hTildePrime/hTilde^2) = -c * h1 * hTildePrime
        x = X[0]
        phase = x - c * t - x0
        sechSqd = old_div(1.0, np.cosh(r * phase)**2)
        hprime = -2.0 * (h2 - h1) * r * sechSqd * np.tanh(r * phase)
        hw = -c * h1 * hprime
        return hw


###############################
#     BOUNDARY CONDITIONS     #
###############################


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
analytical_Solution = {'h_exact': water_height_at_tfinal(),
                       'hu_exact': Zero(),
                       'hv_exact': Zero(),
                       'heta_exact': Zero(),
                       'hw_exact': Zero()}

mySWFlowProblem = SWFlowProblem.SWFlowProblem(sw_model=opts.sw_model,
                                              cfl=0.25,
                                              outputStepping=outputStepping,
                                              structured=True,
                                              he=he,
                                              nnx=nnx,
                                              nny=nny,
                                              domain=domain,
                                              initialConditions=initialConditions,
                                              boundaryConditions=boundaryConditions,
                                              bathymetry=bathymetry_function,
                                              analyticalSolution=analytical_Solution)
mySWFlowProblem.physical_parameters['LINEAR_FRICTION'] = 0
mySWFlowProblem.physical_parameters['mannings'] = 0.0
