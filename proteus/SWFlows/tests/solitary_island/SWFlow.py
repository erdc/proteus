from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2DCV
from proteus.mprans import GN_SW2DCV
from proteus.Domain import RectangularDomain
import numpy as np
from proteus.Gauges import PointGauges
from proteus import (Domain, Context,
                     MeshTools as mt)
from proteus.Profiling import logEvent
import proteus.SWFlows.SWFlowProblem as SWFlowProblem

"""
This is the problem of a solitary wave overtopping a canonical island.
This experiment was done at (insert ref here). There were several gauges
places in the tank measuring the water height (location can be see in bottom)
"""

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts = Context.Options([
    ('sw_model', 1, "sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time", 10.0, "Final time for simulation"),
    ("dt_output", 0.1, "Time interval to output solution"),
    ("cfl", 0.25, "Desired CFL restriction"),
    ("refinement", 4, "Refinement level"),
    ("reflecting_BCs",False,"Use reflecting BCs")
])

###################
# DOMAIN AND MESH #
###################
L = (30.0, 25.0)  # this is length in x direction and y direction
refinement = opts.refinement
domain = RectangularDomain(L=L)

# CREATE REFINEMENT #
nnx0 = 6
nnx = (nnx0 - 1) * (2**refinement) + 1
nny = old_div((nnx - 1), 2) + 1
he = old_div(L[0], float(nnx - 1))
triangleOptions="pAq30Dena%f"  % (0.5*he**2,)

###############################
#  CONSTANTS NEEDED FOR SETUP #
###############################
g = 9.81

# stuff for solitary wave
h0 = 0.32
alpha = 0.181 * h0
xs = 7.56
r = np.sqrt(old_div(3. * alpha, (4. * h0**2 * (h0 + alpha))))
c = np.sqrt(g * (h0 + alpha))

# stuff for cone bathymetry
htop = 0.625
rcone = 3.6
scone = 4.0
hcone = 0.9


#####################################
#   Some functions defined here    #
####################################
def solitary_wave(x, t):
    sechSqd = (1.00 / np.cosh(r * (x - xs)))**2.00
    return alpha * sechSqd


def bathymetry_function(X):
    x = X[0]
    y = X[1]
    radius = np.sqrt((x - 12.96)**2 + (y - 13.80)**2)

    # need to do this annoying thing for piecewise functions
    conds = [radius <= rcone, radius > rcone]
    bath = [lambda radius: np.minimum(htop, hcone-radius/scone), lambda radius: 0.0]
    return np.piecewise(radius, conds, bath)

##############################
##### INITIAL CONDITIONS #####
##############################
class water_height_at_t0(object):
    def uOfXT(self, X, t):
        hTilde = h0 + solitary_wave(X[0], 0)
        h = max(hTilde - bathymetry_function(X), 0.)
        return h


class x_mom_at_t0(object):
    def uOfXT(self, X, t):
        hTilde = h0 + solitary_wave(X[0], 0)
        h = max(hTilde - bathymetry_function(X), 0.)
        return h * c * old_div(hTilde-h0, hTilde)


class y_mom_at_t0(object):
    def uOfXT(self, X, t):
        return 0.


class heta_at_t0(object):
    def uOfXT(self, X, t):
        h = water_height_at_t0().uOfXT(X, t)
        return h**2


class hw_at_t0(object):
    def uOfXT(self, X, t):
        sechSqd = (1.0 / np.cosh(r * (X[0] - xs)))**2.0
        hTilde = h0 + solitary_wave(X[0], 0)
        h = max(hTilde - bathymetry_function(X), 0.)
        hTildePrime = -2.0 * alpha * r * np.tanh(r*(X[0]-xs)) * sechSqd
        hw = -h**2 * old_div(c * h0 * hTildePrime, hTilde**2)
        return hw

###############################
##### BOUNDARY CONDITIONS #####
###############################
X_coords = (0.0, 30.0)  # this is x domain, used in BCs
Y_coords = (0.0, 25.0)  # this is y domain, used in BCs

def x_mom_DBC(X, flag):
    if X[0] == X_coords[0] or X[0] == X_coords[1]:
        return lambda X, t: 0.0


def y_mom_DBC(X, flag):
    if X[1] == Y_coords[0] or X[1] == Y_coords[1]:
        return lambda X, t: 0.0

# **************************** #
# ********** GAUGES ********** #
# **************************** #
want_gauges = False
heightPointGauges = PointGauges(gauges=((('h'), ((7.56, 16.05,  0),
                                 (7.56, 14.55, 0),
                                 (7.56, 13.05, 0),
                                 (7.56, 11.55, 0),
                                 (9.36, 13.80, 0),
                                 (10.36, 13.80, 0),
                                 (12.96, 11.22, 0),
                                 (15.56, 13.80, 0))),),
                activeTime=(0., opts.final_time),
                fileName='wave_gauges.csv')

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
                      'y_mom': y_mom_DBC,
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
mySWFlowProblem.physical_parameters['mannings'] = 0.0
if want_gauges:
    mySWFlowProblem.auxiliaryVariables = [heightPointGauges]
