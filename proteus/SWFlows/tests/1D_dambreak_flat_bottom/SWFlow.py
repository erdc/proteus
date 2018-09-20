from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2D
from proteus.mprans import SW2DCV
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
    ('start_at_t0',False,"Consider dambreak from t=0?"),
    ('refinement',3,"Refinement level"),
    ('sw_model',0,"sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time",5.0,"Final time for simulation"),
    ("dt_output",0.1,"Time interval to output solution"),
    ("cfl",0.33,"Desired CFL restriction")
    ])

###################
# DOMAIN AND MESH #
###################
L=(10.0,1.0)
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
bathymetry = lambda x: 0.*x[0]

##############################
##### INITIAL CONDITIONS #####
##############################
hl=0.005
xc=5.0
g=9.81
class dam_break_problem_starting_at_t0(object):
    def __init__(self):
        pass
    def uOfXT(self,X,t):
        import math
        x = X[0]
        if (x <= xc):
            h = hl
        else:
            h = 0.000
        return h

class dam_break_problem_starting_at_t1(object):
    def __init__(self):
        pass
    def uOfXT(self,X,t):
        import math
        x = X[0]
        xA = xc-math.sqrt(g*hl)
        xB = xc+2*math.sqrt(g*hl)
        if (0 <= x and x <= xA):
            return hl
        elif (xA < x <= xB):
            return 4/9./g*(math.sqrt(g*hl)-old_div((x-xc),2.))**2
        else: 
            return 0.

class velX_starting_at_t1(object):
    def __init__(self):
        pass 
    def uOfXT(self,X,t):
        import math
        x = X[0]
        xA = xc-math.sqrt(g*hl)
        xB = xc+2*math.sqrt(g*hl)
        if (0 <= x and x <= xA):
            return 0.
        elif (xA < x <= xB):
            return 2/3.*(x-xc+math.sqrt(g*hl))
        else: 
            return 0.

class momX_starting_at_t1(object):
    def __init__(self):
        pass
    def uOfXT(self,X,t):
        h = dam_break_problem_starting_at_t1()
        vel = velX_starting_at_t1()
        return h.uOfXT(X,t)*vel.uOfXT(X,t)

class Zero(object):
    def uOfXT(self,x,t):
        return 0.0

# ********************************** #
# ***** Create mySWFlowProblem ***** #
# ********************************** #
outputStepping = SWFlowProblem.OutputStepping(opts.final_time,dt_output=opts.dt_output)
if opts.start_at_t0:
    initialConditions = {'water_height': dam_break_problem_starting_at_t0(),
                         'x_mom': Zero(),
                         'y_mom': Zero()}
else:
    initialConditions = {'water_height': dam_break_problem_starting_at_t1(),
                         'x_mom': momX_starting_at_t1(),
                         'y_mom': Zero()}
#
boundaryConditions = {'water_height': lambda x,flag: None,
                      'x_mom': lambda x,flag: None,
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
