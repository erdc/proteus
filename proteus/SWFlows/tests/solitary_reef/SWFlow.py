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

"""
This is the problem of a solitary wave overtopping a canonical island.
This experiment was done at (insert ref here). There were several gauges
places in the tank measuring the water height (location can be see in bottom)
"""

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts = Context.Options([
    ('sw_model', 0, "sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time", 9.0, "Final time for simulation"),
    ("dt_output", 0.1, "Time interval to output solution"),
    ("cfl", 0.30, "Desired CFL restriction"),
    ("refinement", 4, "Refinement level"),
    ("reflecting_BCs",True,"Use reflecting BCs")
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
h0 = 0.78
alpha = 0.5 * h0
xs = 5.0
vel = np.sqrt(old_div(3 * alpha, (4 * h0**3)))

# stuff for bathymetry, including shelf and cone
rcone = 3.
hcone = 0.45
yc = 13.25


#####################################
#   Some functions defined here    #
####################################
def solitary_wave(x, t):
    sechSqd = (1.00 / np.cosh(vel * (x - xs)))**2.00
    soliton = alpha * sechSqd
    return soliton


def bathymetry_function(X):
    x = X[0]
    y = X[1]
    # first define cone topography
    cone = np.max(hcone - np.sqrt(((x-17.0)**2+(y-yc)**2)/(rcone/hcone)**2),0.0)

    # define some stuff we need for shelf
    dist = 1.0 - MIN(1.0, ABS(y - yc)/yc)
    aux_x = 12.50 + 12.4999 * (1.0 - dist)
    aux_z = 0.70 + 0.050 * (1.0 - dist)

    # then do annoying stuff for piecewise functions
    conds = [x < 10.2, (10.2 < x) & (x - aux_x <= 0.0), (aux_x - x <= 0.) & (x <= 25.), \
                        (13.948 < x) & (x<= 14.045)]
    bath = [lambda x: 0, \
    lambda x: aux_z / (aux_x-10.20) * (x-10.2), \
    lambda x: 0.75 + (aux_z-0.75)/(aux_x-25.)*(x-25.), \
    lambda x: 3.6 / 20. + 0.076 - (0.076 - 0.022) / (14.045 - 13.948) * (x - 13.948) - h0, \
    lambda x: 1 / 20. * (x - 10.) - h0]
    return np.piecewise(x, conds, bath)

##############################
##### INITIAL CONDITIONS #####
##############################
class water_height_at_t0(object):
    def uOfXT(self, X, t):
        eta = h0 + solitary_wave(X[0], 0)
        h = max(eta - bathymetry_function(X), 0.)
        return h


class x_mom_at_t0(object):
    def uOfXT(self, X, t):
        return solitary_wave(X[0], 0) * np.sqrt(g / h0)


class y_mom_at_t0(object):
    def uOfXT(self, X, t):
        return 0.


class heta_at_t0(object):
    def uOfXT(self, X, t):
        h = water_height_at_t0().uOfXT(X, t)
        return h**2


class hw_at_t0(object):
    def uOfXT(self, X, t):
        hw = 2.0 * solitary_wave(X[0], 0) * vel * np.sinh(vel * (X[0]-xs))
        return hw

###############################
##### BOUNDARY CONDITIONS #####
###############################

# Actually don't need any of these
def water_height_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda x, t: water_height_at_t0().uOfXT(X, 0.0)

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

# **************************** #
# ********** GAUGES ********** #
# **************************** #
from proteus.Gauges import PointGauges
p = PointGauges(gauges=(( ('h'), ((30, 15, 0), (60,15,0)) ),),
                activeTime=(0, 10),
                sampleRate=0.1,
                fileName='island_gauges.csv')
auxiliaryVariables = [p]


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
boundaryConditions = {'water_height': lambda x,flag: None,
                      'x_mom': lambda x,flag: None,
                      'y_mom': lambda x,flag: None,
                      'h_times_eta': lambda x,flag: None,
                      'h_times_w': lambda x,flag: None}
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
