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


# Note that this set up for a solitary wave is from the paper
# A rapid numerical method for solving Serre Green Naghdi
# equations describing long free surface gravity waves by
# Favrie and Gavrilyuk equation 31. It assumes the reference height is h1.

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
L = (20.0, 1.0)
refinement = opts.refinement
domain = RectangularDomain(L=L, x=[0, 0, 0])

# CREATE REFINEMENT #
nnx0 = 6
nnx = (nnx0 - 1) * (2**refinement) + 1
nny = old_div((nnx - 1), 10) + 1

he = old_div(L[0], float(nnx - 1))
triangleOptions = "pAq30Dena%f" % (0.5 * he**2,)

#################
# SOLITARY WAVE #
#################
h1 = 0.1
h2 = 0.11
x0 = 2.0 # initial location of solitary wave
g = 9.81


def soliton(x, t):
    D = np.sqrt(g * h2)
    z = np.sqrt(old_div(3.0 * (h2 - h1), h2 * h1**2))
    phase = x - D * t - x0
    a1 = z * phase / 2.0
    return h1 + (h2 - h1) * 1.0 / np.cosh(a1)**2


def u(x, t):
    D = np.sqrt(g * h2)
    return D * (1.0 - old_div(h1, soliton(x, t)))

###############################
##### BOUNDARY CONDITIONS #####
###############################


def water_height_DBC(X, flag):
    if X[0] == 0:
        return lambda x, t: soliton(X[0], t)
    elif X[0] == L[0]:
        return lambda X, t: h1


def x_mom_DBC(X, flag):
    if X[0] == 0:
        return lambda X, t: soliton(X[0], t) * u(X[0], t)
    elif X[0] == L[0]:
        return lambda X, t: 0.0

##############################
##### INITIAL CONDITIONS #####
##############################
class water_height_at_t0(object):
    def uOfXT(self, X, t):
        return soliton(X[0], t)

class x_mom_at_t0(object):
    def uOfXT(self, X, t):
        return soliton(X[0], t) * u(X[0],t)

class Zero(object):
    def uOfXT(self, X, t):
        return 0.0


# ********************************** #
# ***** Create mySWFlowProblem ***** #
# ********************************** #
outputStepping = SWFlowProblem.OutputStepping(
    opts.final_time, dt_output=opts.dt_output)
initialConditions = {'water_height': water_height_at_t0(),
                     'x_mom': x_mom_at_t0(),
                     'y_mom': Zero()}
boundaryConditions = {'water_height': water_height_DBC,
                      'x_mom': x_mom_DBC,
                      'y_mom': lambda x, flag: lambda x, t: 0.0}
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
                                              bathymetry=lambda X: X[0] * 0)
mySWFlowProblem.physical_parameters['LINEAR_FRICTION'] = 0
mySWFlowProblem.physical_parameters['mannings'] = 0
