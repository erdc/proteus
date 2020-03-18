from __future__ import division
from builtins import object
from past.utils import old_div
from proteus.mprans import (SW2DCV, GN_SW2DCV)
from proteus.Domain import RectangularDomain
import numpy as np
from proteus import (Domain, Context,
                     MeshTools as mt)
from proteus.Profiling import logEvent
import proteus.SWFlow.SWFlowProblem as SWFlowProblem

"""This is a dam-break problem with three obstacles.  The computation domain is
D = [0,75m]x[0,30m]. The dam is located at x = 16m and is of height 1.875m. """

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts = Context.Options([
    ('sw_model', 0, "sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time", 30.0, "Final time for simulation"),
    ("dt_output", 0.1, "Time interval to output solution"),
    ("refinement", 4, "Level of refinement"),
    ("cfl", 0.33, "Desired CFL restriction"),
    ("reflecting_BCs", True, "Use reflecting BCs")
])

###################
# DOMAIN AND MESH #
###################
L = (75.0, 30.0)
refinement = opts.refinement
domain = RectangularDomain(L=L)

# CREATE REFINEMENT #
nnx0 = 6
nnx = (nnx0 - 1) * (2**refinement) + 1
nny = old_div((nnx - 1), 2) + 1
he = old_div(L[0], float(nnx - 1))
triangleOptions = "pAq30Dena%f" % (0.5 * he**2,)

######################
##### BATHYMETRY #####
######################


def bathymetry_function(X):
    x = X[0]
    y = X[1]
    bump1 = 1 - 1. / 8 * np.sqrt((x - 30)**2 + (y - 6)**2)
    bump2 = 1 - 1. / 8 * np.sqrt((x - 30)**2 + (y - 24)**2)
    bump3 = 3 - 3. / 10 * np.sqrt((x - 47.5)**2 + (y - 15)**2)
    return np.maximum(np.maximum(np.maximum(0., bump1), bump2), bump3)

##############################
##### INITIAL CONDITIONS #####
##############################


class water_height_at_t0(object):
    def uOfXT(self, X, t):
        x = X[0]
        if (x <= 16):
            eta = 1.875
        else:
            eta = 0.

        z = bathymetry_function(X)
        return max(eta - z, 0.)


class Zero(object):
    def uOfXT(self, x, t):
        return 0.0

class heta_at_t0(object):
    def uOfXT(self, X, t):
        h = water_height_at_t0().uOfXT(X, t)
        return h**2


class hw_at_t0(object):
    def uOfXT(self, X, t):
        return 0.0

###############################
##### BOUNDARY CONDITIONS #####
###############################
# If we want other BCs instead of reflecting
X_coords = (0.0, L[0])  # this is x domain, used in BCs
Y_coords = (0.0, L[1])  # this is y domain, used in BCs

def x_mom_DBC(X, flag):
    if X[0] == X_coords[0] or X[0] == X_coords[1]:
        return lambda X, t: 0.0


def y_mom_DBC(X, flag):
    if X[1] == Y_coords[0] or X[1] == Y_coords[1]:
        return lambda X, t: 0.0

# ********************************** #
# ***** Create mySWFlowProblem ***** #
# ********************************** #
outputStepping = SWFlowProblem.OutputStepping(
    opts.final_time, dt_output=opts.dt_output)
initialConditions = {'water_height': water_height_at_t0(),
                     'x_mom': Zero(),
                     'y_mom': Zero(),
                     'h_times_eta': heta_at_t0(),
                     'h_times_w': Zero()}
boundaryConditions = {'water_height': lambda x, flag: None,
                      'x_mom': lambda x, flag: None,
                      'y_mom': lambda x, flag: None,
                      'h_times_eta': lambda x, flag: None,
                      'h_times_w': lambda x, flag: None}
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
mySWFlowProblem.physical_parameters['LINEAR_FRICTION'] = 0
mySWFlowProblem.physical_parameters['mannings'] = 0.02
