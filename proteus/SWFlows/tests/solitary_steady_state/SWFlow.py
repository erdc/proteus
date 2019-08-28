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
This is the problem of steady state solutions to the green-naghdi equations.
We use the Bernoulli relation to derive the bathymetry profile where
we assume the water height is given by h = h0 + a h0 sech(r*(x-x0))^2.
The variables are defined below. See 'Robust Explicit Relaxtion Technique
For Solving The Green-Naghdi Equations' by Guermond, Kees, Popov, Tovar.
Note that this is a fake 1D problem, ie
we are doing simulation in 2d but only consider x direction velocity
"""

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
g = 9.81  # define gravity
opts = Context.Options([
    ('sw_model', 1, "sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time", 100, "Final time for simulation"),
    ("dt_output", 0.1, "Time interval to output solution"),
    ("cfl", 0.20, "Desired CFL restriction"),
    ("refinement", 4, "Refinement level")
])

###################
# DOMAIN AND MESH #
###################
L = (20.0, 1.0)
refinement = opts.refinement
domain = RectangularDomain(L=L, x=[-10.0, 0, 0])
X_coords = (-10.0, 10.0)  # this is domain, used in BCs

# CREATE REFINEMENT #
nnx0 = 6
nnx = (nnx0 - 1) * (2**refinement) + 1
nny = old_div((nnx - 1), 10) + 1
he = old_div(L[0], float(nnx - 1))
triangleOptions = "pAq30Dena%f" % (0.5 * he**2,)

###############################
#  CONSTANTS NEEDED FOR SETUP #
###############################
q = 1  # this is flow rate
a = 0.2  # this is amplitude
r = 1  # this is solitary wave width
# we define the reference height in terms of r
h0 = np.sqrt(old_div(3 * a, (4 * r**2*(1+a))))
cBer = h0 + old_div(q**2 , (2 * g * h0**2))  # this is Bernoulli constant
x0 = 0  # wave is centered at x = 0

###############################
#   Functions defined here    #
###############################
def solitary_wave(x, t):
    sechSqd = old_div(1.0, np.cosh(r*(x-x0))**2.0)
    soliton = h0 + a * h0 * sechSqd
    return soliton


def bathymetry_function(X):
    x = X[0]
    num = -3 * 3**(1/3) * a**(3/2) * g + 8 * q**2 * r**2 * np.sqrt((1+a)*r**2)
    denom = 6 * g * np.sqrt((1+a)*r**2)
    sechSqd = old_div(1.0, np.cosh(r*(x-x0))**2.0)
    z = old_div(num, denom) * sechSqd
    return z

##############################
##### INITIAL CONDITIONS #####
##############################
class water_height_at_t0(object):
    def uOfXT(self, X, t):
        h = h0 + a * h0 * solitary_wave(X[0], 0)
        return h

class x_mom_at_t0(object):
    def uOfXT(self, X, t):
        return 0.95 * q

class x_mom_exact(object):
    def uOfXT(self, X, t):
        return q

class y_mom_at_t0(object):
    def uOfXT(self, X, t):
        return 0.

# For analytical solution
class Zero(object):
    def uOfXT(self, x, t):
        return 0.0

"""
heta and hw are needed for the dispersive modified green naghdi equations
source is 'ROBUST EXPLICIT RELAXATION TECHNIQUE FOR SOLVING
THE GREEN NAGHDI EQUATIONS' by Guermond, Kees, Popov, Tovar
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
def water_height_DBC(X, flag):
    if (opts.sw_model==1):
        if X[0] == X_coords[0]:
            return lambda x, t: water_height_at_t0().uOfXT(X, 0.0)
        elif X[0]==X_coords[1]:
            return lambda x,t: water_height_at_t0().uOfXT(X ,0.0)
    if (opts.sw_model==0):
        if X[0] == X_coords[1]:
            return lambda x, t: water_height_at_t0().uOfXT(X, 0.0)

def x_mom_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda X, t: x_mom_exact().uOfXT(X, 0.0)


def y_mom_DBC(X, flag):
    return lambda x, t: 0.0


def heta_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda x, t: heta_at_t0().uOfXT(X, 0.0)
    elif X[0]==X_coords[1]:
        return lambda x,t: heta_at_t0().uOfXT(X, 0.0)


def hw_DBC(X, flag):
    if X[0] == X_coords[0]:
        return lambda x, t: hw_at_t0().uOfXT(X, 0.0)
    elif X[0]==X_coords[1]:
        return lambda x,t: hw_at_t0().uOfXT(X, 0.0)


# ********************************** #
# ***** Create mySWFlowProblem ***** #
# ********************************** #
outputStepping = SWFlowProblem.OutputStepping(opts.final_time, dt_output=opts.dt_output)

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
analyticalSolution={'h_exact': water_height_at_t0(),
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
mySWFlowProblem.physical_parameters['mannings'] = 0
