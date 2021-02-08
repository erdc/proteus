from __future__ import division
from builtins import object
from past.utils import old_div
from proteus.mprans import (SW2DCV, GN_SW2DCV)
from proteus.Domain import RectangularDomain, PlanarStraightLineGraphDomain
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
    ("final_time", 20.0, "Final time for simulation"),
    ("dt_output", 0.1, "Time interval to output solution"),
    ("refinement", 4, "Level of refinement"),
    ("structured", True, "Structured or unstructured mesh"),
    ("he", 0.5, "Mesh size for unstructured mesh"),
    ("cfl", 0.33, "Desired CFL restriction"),
    ("reflecting_BCs", True, "Use reflecting BCs"),
    ("mannings", 0.02, "Mannings roughness coefficient")
])

###################
# DOMAIN AND MESH #
###################
L = (75.0, 30.0)
refinement = opts.refinement
rectangle = RectangularDomain(L=L)

# CREATE REFINEMENT #
nnx0 = 6
nnx = (nnx0 - 1) * (2**refinement) + 1
nny = old_div((nnx - 1), 2) + 1
he = old_div(L[0], float(nnx - 1))
if opts.structured:
    domain = rectangle
else:
    rectangle.writePoly("dam3bumps")
    domain = PlanarStraightLineGraphDomain(fileprefix="dam3bumps")
    domain.MeshOptions.triangleOptions = "pAq30Dena%f" % (0.5 * opts.he**2,)
    nnx = None
    nny = None

######################
##### BATHYMETRY #####
######################


def bathymetry_function(X):
    x = X[0]
    y = X[1]
    bump1 = 1. - 1. / 8. * np.sqrt((x - 30.)**2 + (y - 6.)**2)
    bump2 = 1. - 1. / 8. * np.sqrt((x - 30.)**2 + (y - 24.)**2)
    bump3 = 3. - 3. / 10. * np.sqrt((x - 47.5)**2 + (y - 15.)**2)
    return np.maximum(np.maximum(np.maximum(0., bump1), bump2), bump3)


##############################
##### INITIAL CONDITIONS #####
##############################

class water_height_at_t0(object):
    def uOfXT(self, X, t):
        x = X[0]
        if (x <= 16.0):
            eta = 1.875
        else:
            eta = 0.0
        return max(eta, 0.)


class Zero(object):
    def uOfXT(self, x, t):
        return 0.0

"""
heta and hw are needed for the hyperbolic serre-green-naghdi equations.
For initial conditions, heta -> h^2, hbeta->q(dot)grad(Z), hw -> h^2div(u)+3/2*hbeta.
It's often okay to take hbeta=0. Note that the BCs for the heta and hw should be same as h
and BCs for hbeta should be same as x_mom.
For more details see: 'Hyperbolic relaxation technique for solving the dispersive Serre Equations
with topography' by Guermond, Popov, Tovar, Kees.
"""        

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
                     'h_times_w': Zero(),
                     'h_times_beta': Zero()}
boundaryConditions = {'water_height': lambda x, flag: None,
                      'x_mom': x_mom_DBC,
                      'y_mom': y_mom_DBC,
                      'h_times_eta': lambda x, flag: None,
                      'h_times_w': lambda x, flag: None,
                      'h_times_beta': x_mom_DBC}
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
                                              reflectingBCs=opts.reflecting_BCs,
                                              bathymetry=bathymetry_function)
mySWFlowProblem.physical_parameters['LINEAR_FRICTION'] = 0
mySWFlowProblem.physical_parameters['mannings'] = opts.mannings
mySWFlowProblem.physical_parameters['cE'] = 1.0
