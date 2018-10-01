from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.Domain import RectangularDomain
import numpy as np
from proteus import (Domain, Context,
                     MeshTools as mt)
from proteus.Profiling import logEvent
import proteus.SWFlows.SWFlowProblem as SWFlowProblem

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ('sw_model',0,"sw_model = {0,1} for {SWEs,DSWEs}"),
    ('refinement',3,"Refinement level"),
    ("final_time",1000.0,"Final time for simulation"),
    ("dt_output",100.0,"Time interval to output solution"),
    ("cfl",0.33,"Desired CFL restriction")
    ])

###################
# DOMAIN AND MESH #
###################
L=(8000.0,800.0)
refinement = opts.refinement
domain = RectangularDomain(L=L,x=[0,0,0])

# CREATE REFINEMENT #
nnx0=6
nnx = (nnx0-1)*(2**refinement)+1
nny = old_div((nnx-1),10)+1

he = old_div(L[0],float(nnx-1))
triangleOptions="pAq30Dena%f"  % (0.5*he**2,)

######################
##### BATHYMETRY #####
######################
h0=10
a=3000
B=2
k=0.001
g = SWFlowProblem.default_physical_parameters['gravity']
p = old_div(np.sqrt(8*g*h0),a)
s = old_div(np.sqrt(p**2 - k**2),2.)
mannings=k

def bathymetry_function(X):
    x=X[0]
    return h0*(x-old_div(L[0],2))**2/a/a

def eta_function(x,t):
    coeff1 = a**2*B**2/8./g/g/h0
    coeff2 = -B**2/4./g
    coeff3 = old_div(-1.,g)

    eta_part1 = coeff1*np.exp(-k*t)*(-s*k*np.sin(2*s*t)+(old_div(k**2,4.)-s**2)*np.cos(2*s*t))
    eta_part2 = coeff2*np.exp(-k*t)
    eta_part3 = coeff3*np.exp(-k*t/2.)*(B*s*np.cos(s*t)+k*B/2.*np.sin(s*t))*(x-old_div(L[0],2))

    return h0 + eta_part1 + eta_part2 + eta_part3

##############################
##### INITIAL CONDITIONS #####
##############################
class water_height_at_t0(object):
    def uOfXT(self,X,t):
        eta = eta_function(X[0],0)
        bathymetry = bathymetry_function(X)
        h = eta-bathymetry
        hp = max(h,0.)
        return hp

class Zero(object):
    def uOfXT(self,X,t):
        return 0.0

##########################
##### EXACT SOLUTION #####
##########################
class water_height(object):
    def uOfXT(self,X,t):
        h = eta_function(X[0],t) - bathymetry_function(X)
        hp = max(h,0.)
        return hp

class exact_momx(object):
    def uOfXT(self,X,t):
        ht = eta_function(X[0],t) - bathymetry_function(X)
        if (ht >= 0):
            h=ht
        else:
            h=0
        u = B*np.exp(-k*t/2.)*np.sin(s*t)
        return h*u

analyticalSolution = {0: water_height(),
                      1: exact_momx(),
                      2: Zero()}

# ********************************** #
# ***** Create mySWFlowProblem ***** #
# ********************************** #
outputStepping = SWFlowProblem.OutputStepping(opts.final_time,
                                              dt_output=opts.dt_output,
                                              dt_init=1E-6)
initialConditions = {'water_height': water_height_at_t0(),
                     'x_mom': Zero(),
                     'y_mom': Zero()}
boundaryConditions = {'water_height': lambda x,flag: None,
                      'x_mom': lambda x,flag: None,
                      'y_mom': lambda x,flag: lambda x,t: 0.0}
mySWFlowProblem = SWFlowProblem.SWFlowProblem(sw_model=0,
                                              cfl=opts.cfl,
                                              outputStepping=outputStepping,
                                              structured=True,
                                              triangleFlag=0,
                                              he=he,
                                              nnx=nnx,
                                              nny=nny,
                                              domain=domain,
                                              initialConditions=initialConditions,
                                              boundaryConditions=boundaryConditions,
                                              bathymetry=bathymetry_function,
                                              analyticalSolution=analyticalSolution)
mySWFlowProblem.physical_parameters['LINEAR_FRICTION']=1
mySWFlowProblem.physical_parameters['mannings']=mannings
