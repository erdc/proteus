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

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ('sw_model',0,"sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time",25.0,"Final time for simulation"),
    ("dt_output",0.1,"Time interval to output solution"),
    ("cfl",0.33,"Desired CFL restriction"),
    ("refinement",4,"Refinement level")
    ])

###################
# DOMAIN AND MESH #
###################
L=(100.0,10.) #45,4.5
L0=[-90,0,0]  #-35,0,0
refinement = opts.refinement
domain = RectangularDomain(L=L,x=L0)

# CREATE REFINEMENT #
nnx0=6
nnx = (nnx0-1)*(2**refinement)+1
nny = old_div((nnx-1),10)+1

he = old_div(L[0],float(nnx-1))
triangleOptions="pAq30Dena%f"  % (0.5*he**2,)

#################
# SOLITARY WAVE #
#################
g = SWFlowProblem.default_physical_parameters['gravity']
h0 = 1.0
a = 0.30  # amplitude
slope = 1.0 / 19.850
k_wavenumber = np.sqrt(3.0 * a/(4.0 * h0**3))  # wavenumber
z = np.sqrt(3.0 * a * h0) / (2.0 * h0 * np.sqrt(h0 * (1.0 + a)))
wavelength = 2.0 / k_wavenumber * np.arccosh(np.sqrt(1.0 / 0.050))  # wavelength of solitary wave
c = np.sqrt(g * (1.0 + a) * h0)
x0 = -50.0 - h0/slope - wavelength/2.0  # location of the toe of the beach

def soliton(x,t): #
    sechSqd = (1.00/np.cosh( z*(x-x0-c*t)))**2.00
    return a * h0 * sechSqd

def u(x,t):
    eta = soliton(X[0],t)
    return c * eta / (h0 + eta)

def bathymetry(X):
    x=X[0]
    return numpy.maximum(slope*x,-h0)
                                                                           
###############################
##### BOUNDARY CONDITIONS #####
###############################
def water_height_DBC(X,flag):
    if X[0]==L0[0]:
        return lambda x,t: h0
    
def x_mom_DBC(X,flag):
    if X[0]==L0[0]:
        return lambda X,t: 0.0
    
##############################
##### INITIAL CONDITIONS #####
##############################
class water_height_at_t0(object):
    def uOfXT(self,X,t):
        eta = soliton(X[0],0)
        h = eta-bathymetry(X)
        hp = max(h,0.)
        return hp

class x_mom_at_t0(object):
    def uOfXT(self,X,t):
        eta = soliton(X[0],0)
        h = eta-bathymetry(X)
        hp = max(h,0.)
        Umom = hp * c * eta / (h0 + eta)
        return Umom
    
class Zero():
    def uOfXT(self,X,t):
        return 0.0
    
# ********************************** #
# ***** Create mySWFlowProblem ***** #
# ********************************** #
outputStepping = SWFlowProblem.OutputStepping(opts.final_time,dt_output=opts.dt_output)
initialConditions = {'water_height': water_height_at_t0(),
                     'x_mom': x_mom_at_t0(),
                     'y_mom': Zero()}
boundaryConditions = {'water_height': water_height_DBC,
                      'x_mom': x_mom_DBC,
                      'y_mom': lambda x,flag: lambda x,t: 0.0}
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
                                              bathymetry=bathymetry)
mySWFlowProblem.physical_parameters['LINEAR_FRICTION']=0
mySWFlowProblem.physical_parameters['mannings']=0
