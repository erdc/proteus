from __future__ import division
from builtins import object
from past.utils import old_div
from proteus.mprans import (SW2DCV, GN_SW2DCV)
from proteus.Domain import RectangularDomain, PlanarStraightLineGraphDomain
import numpy as np
from proteus import (Domain, Context, MeshTools as mt)
from proteus.Profiling import logEvent
import proteus.SWFlow.SWFlowProblem as SWFlowProblem


"""
This is a simple benchmark of a solitary wave propagating over a flat bottom.
"""

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #

opts = Context.Options([
    ('sw_model', 0, "sw_model = {0,1} for {SWEs, Disperisve SWEs}}"),
    ("final_time", 40.0, "Final time for simulation"),
    ("dt_output", 0.1, "Time interval to output solution"),
    ("cfl", 0.25, "Desired CFL restriction"),
    ("refinement", 4, "Refinement level"),
    ("reflecting_BCs", False, "Use reflecting BCs for all boundaries"),
    ("structured", True, "Structured or unstructured mesh"),
    ("he", 0.1, "Mesh size for unstructured mesh")
])

###################
# DOMAIN AND MESH #
###################
L = (25.0, 1.0)
refinement = opts.refinement
rectangle = RectangularDomain(L=L, x=[0, 0, 0])

# CREATE REFINEMENT #
nnx0 = 6
nnx = (nnx0 - 1) * (2**refinement) + 1
nny = old_div((nnx - 1), 10) + 1
he = old_div(L[0], float(nnx - 1))
triangleOptions = "pAq30Dena%f" % (0.5 * he**2,)
if opts.structured:
    domain = rectangle
else:
    rectangle.writePoly("bump")
    domain = PlanarStraightLineGraphDomain(fileprefix="bump")
    domain.MeshOptions.triangleOptions = "pAq30Dena%f" % (0.5 * opts.he**2,)
    nnx = None
    nny = None

domain.MeshOptions.nnx = nnx
domain.MeshOptions.nny = nny


##########################################
# DEFINE INITIAL CONSTANTS AND FUNCTIONS #
##########################################
g = 9.81

# constants for transcritical bump
hL = 0.28205279813802181
q_in = 0.18

def bathymetry_function(X):
    x = X[0]
    #  define conditionals for bathymetry
    conds = [(8.0 < x) & (x <= 12.)]
    #  define the functions
    bath = [lambda x: 0.2/64.0 * (x - 8.0)**3 * (12.0 - x)**3]
    return np.piecewise(x, conds, bath)


##############################
#    INITIAL CONDITIONS      #
##############################


class water_height_at_t0(object):
    def uOfXT(self, X, t):
        h = max(hL - bathymetry_function(X), 0.)
        return h


class x_mom_at_t0(object):
    def uOfXT(self, X, t):
        return q_in


class y_mom_at_t0(object):
    def uOfXT(self, X, t):
        return 0.


"""
heta and hw are needed for the modified green naghdi equations.
For initial conditions, heta -> h^2 and hw -> h^2*div(u).
Note that the BC flags for the heta and hw should be same as h.
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
        return 0.

class Zero(object):
    def uOfXT(self, x, t):
        return 0.0

###############################
#     BOUNDARY CONDITIONS     #
###############################
X_coords = (0.0, L[0])
Y_coords = (0.0, L[1])

def h_DBC(X, flag):
    if X[0] == X_coords[1]:
        return lambda x, t: hL

def x_mom_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda x, t: q_in

def heta_DBC(X, flag):
    if X[0] == X_coords[1]:
        return lambda x, t: hL**2

def hw_DBC(X, flag):
    if X[0] == X_coords[1]:
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
boundaryConditions = {'water_height': h_DBC,
                      'x_mom': x_mom_DBC,
                      'y_mom': lambda x, flag: lambda x, t: 0.0,
                      'h_times_eta': heta_DBC,
                      'h_times_w': hw_DBC}

mySWFlowProblem = SWFlowProblem.SWFlowProblem(sw_model=opts.sw_model,
                                              cfl=opts.cfl,
                                              outputStepping=outputStepping,
                                              structured=opts.structured,
                                              he=he,
                                              nnx=nnx,
                                              nny=nny,
                                              domain=domain,
                                              initialConditions=initialConditions,
                                              boundaryConditions=boundaryConditions,
                                              bathymetry=bathymetry_function)
mySWFlowProblem.physical_parameters['LINEAR_FRICTION'] = 0
mySWFlowProblem.physical_parameters['mannings'] = 0.0
