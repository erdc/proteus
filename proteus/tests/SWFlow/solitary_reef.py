from __future__ import division
from builtins import object
from past.utils import old_div
from proteus.mprans import (SW2DCV, GN_SW2DCV)
from proteus.Domain import RectangularDomain
import numpy as np
from proteus import (Domain, Context,
                     MeshTools as mt)
from proteus.Profiling import logEvent
from proteus.Gauges import PointGauges
import proteus.SWFlow.SWFlowProblem as SWFlowProblem

"""
We reproduce the 2009-2010 experiments of [Swigler, 2009] and
[Lynett, 2019] performed at the O.H. Hinsdale Wave Research
Laboratory of Oregon State University. The experiments were conducted
to study specific phenomena that are known to occur when solitary
waves propagate over irregular bathymetry such as shoaling,
refraction, breaking, etc. In the experiment, nine wave gauges (WGs)
were placed along the basin to capture the free surface elevation
along with three Acoustic Doppler Velocimeters (ADVs) that
measured the velocities in both horizontal directions.
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
    ("reflecting_BCs", False, "Use reflecting BCs")
])

###################
# DOMAIN AND MESH #
###################
L = (48.8, 26.5)  # this is length in x direction and y direction
refinement = opts.refinement
domain = RectangularDomain(L=L, x=[0, -13.25, 0])

# CREATE REFINEMENT #
nnx0 = 6
nnx = (nnx0 - 1) * (2**refinement) + 1
nny = old_div((nnx - 1), 2) + 1
he = old_div(L[0], float(nnx - 1))
triangleOptions = "pAq30Dena%f" % (0.5 * he**2,)

###############################
#  CONSTANTS NEEDED FOR SETUP #
###############################
g = 9.81

# stuff for solitary wave
h0 = 0.78
alpha = 0.4  # 0.5 * h0
xs = 5.0
r = np.sqrt(old_div(3. * alpha, (4. * h0**2 * (h0 + alpha))))
c = np.sqrt(g * (h0 + alpha))

# stuff for bathymetry, including shelf and cone
rcone = 3.
hcone = 0.45
yc = 13.25


#####################################
#   Some functions defined here    #
####################################

def solitary_wave(x, t):
    sechSqd = (1.0 / np.cosh(r * (x - xs - c * t)))**2
    return alpha * sechSqd


def bathymetry_function(X):
    x = X[0]
    y = X[1] + yc

    # first define cone topography
    cone = np.maximum(
        hcone - np.sqrt(((x - 17.0)**2 + (y - yc)**2) / (rcone / hcone)**2), 0.0)

    # define piecewise function for base
    base = 0. * x
    conds = [x < 10.2, (10.2 < x) & (x <= 17.5), (17.5 <= x) & (x <= 32.5),
                        32.5 < x]
    base_values = [lambda x: 0.0,
                   lambda x: (0.5 - 0.0) / (17.5 - 10.20) * (x - 10.2),
                   lambda x:  1.0 + (1.0 - 0.5)/(32.5 - 17.5) * (x - 32.5),
                   lambda x: 1.]

    base = np.piecewise(x, conds, base_values)

    # define  piecewise function for shelf
    shelf = 0. * x
    dist = 1.0 - np.minimum(1.0, np.abs(y - yc) / yc)
    aux_x = 12.50 + 12.4999 * (1.0 - dist)
    aux_z = 0.70 + 0.050 * (1.0 - dist)

    conds = [x < 10.2, (10.2 <= x) & (x <= aux_x), (aux_x <= x) & (x <= 25.),
            (25. < x) & (x <= 32.5), 32.5 < x]
    shelf_values = [0.0,
                    aux_z / (aux_x - 10.20) * (x - 10.2),
                    0.75 + (aux_z - 0.75) / (aux_x - 25.) * (x - 25.),
                    1. + (1. - 0.5) / (32.5 - 17.5) * (x - 32.5),
                    1.]
    shelf = np.select(conds, shelf_values)

    bath = np.maximum(base, shelf) + cone
    return bath


######################
# INITIAL CONDITIONS #
######################


class water_height_at_t0(object):
    def uOfXT(self, X, t):
        hTilde = h0 + solitary_wave(X[0], 0)
        h = max(hTilde - bathymetry_function(X), 0.)
        return h


class x_mom_at_t0(object):
    def uOfXT(self, X, t):
        hTilde = h0 + solitary_wave(X[0], 0)
        h = max(hTilde - bathymetry_function(X), 0.)
        return h * c * old_div(hTilde - h0, hTilde)


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
        hTildePrime = -2.0 * alpha * r * np.tanh(r * (X[0] - xs)) * sechSqd
        hw = -h**2 * old_div(c * h0 * hTildePrime, hTilde**2)
        return hw


###############################
##### BOUNDARY CONDITIONS #####
###############################
X_coords = (0.0, 48.8)  # this is x domain, used in BCs
Y_coords = (-13.25, 13.25)  # this is y domain, used in BCs


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
                     'x_mom': x_mom_at_t0(),
                     'y_mom': y_mom_at_t0(),
                     'h_times_eta': heta_at_t0(),
                     'h_times_w': hw_at_t0()}
boundaryConditions = {'water_height': lambda x, flag: None,
                      'x_mom': x_mom_DBC,
                      'y_mom': y_mom_DBC,
                      'h_times_eta': lambda x, flag: None,
                      'h_times_w': lambda x, flag: None}
# **************************** #
# ********** GAUGES ********** #
# **************************** #
want_gauges = False
heightPointGauges = PointGauges(gauges=((('h'), ((7.5, 0.0,  0),
                                                 (13.0, 0.0, 0),
                                                 (21.0, 0.0, 0),
                                                 (7.5, 5.0, 0),
                                                 (13.0, 5.0, 0),
                                                 (21.0, 5.0, 0),
                                                 (25.0, 0.0, 0),
                                                 (25.0, 5.0, 0),
                                                 (25.0, 10.0, 0))),),
                                activeTime=(0., opts.final_time),
                                fileName='reef_wave_gauges.csv')

# ********************************************* #
# ********** Create my SWFlowProblem ********** #
# ********************************************* #
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
                                              bathymetry=bathymetry_function,
                                              analyticalSolution=None)
mySWFlowProblem.physical_parameters['LINEAR_FRICTION'] = 0
mySWFlowProblem.physical_parameters['mannings'] = 0.0
if want_gauges:
    mySWFlowProblem.auxiliaryVariables = [heightPointGauges]
