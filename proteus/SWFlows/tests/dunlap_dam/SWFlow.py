from __future__ import division
from builtins import object
from past.utils import old_div
from proteus import *
from proteus.default_p import *
from proteus.mprans import SW2D
from proteus.mprans import SW2DCV
from proteus.mprans import GN_SW2DCV
from proteus.Domain import RectangularDomain
from proteus.Domain import Mesh2DMDomain
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
    ("final_time",2.0,"Final time for simulation"),
    ("dt_output",1.0,"Time interval to output solution"),
    ("cfl",0.33,"Desired CFL restriction"),
    ("refinement",4,"Refinement level"),
    ("reflectingBCs",True,"Use reflecting BCs")
    ])

###################
# DOMAIN AND MESH #
###################
AdH_file = "Lake_Dunlap_UTM_15N_Meters"
domain = None

##############
# BATHYMETRY #
##############
bathymetry = None
#using mesh z coord

###############################
##### BOUNDARY CONDITIONS #####
###############################
# REFLECTING BCs

##############################
##### INITIAL CONDITIONS #####
##############################
class water_height_at_t0(object):
    """set the water level to 100m behind the dam and dry elsewhere"""
    def uOfXT(self,X,t):
        x = X[0]
        #LINE 1
        x1 = 4701.18
        y1 = 4143.41
        x2 = 4655.5
        y2 = 4392.1
        m = old_div((y2-y1),(x2-x1))
        dam1 = m*(x-x1)+y1

        #LINE 2
        x1 = 4655.5
        y1 = 4392.1
        x2 = 4000.0
        y2 = 5500.0
        m = old_div((y2-y1),(x2-x1))
        dam2 = m*(x-x1)+y1

        return 1.

        # if (X[1] <= dam1 and X[1] <= dam2):
        #     return np.maximum(100.0-X[2],0.)
        # else:
        #     return 0.

class Zero(object):
    """still water conditions"""
    def uOfXT(self,x,t):
        return 0.0

# ********************************** #
# ***** Create mySWFlowProblem ***** #
# ********************************** #
outputStepping = SWFlowProblem.OutputStepping(opts.final_time,dt_output=opts.dt_output)
initialConditions = {'water_height': water_height_at_t0(),
                     'x_mom': Zero(),
                     'y_mom': Zero()}
boundaryConditions = {'water_height': lambda x,flag: None,
                      'x_mom': lambda x,flag: None,
                      'y_mom': lambda x,flag: None}
mySWFlowProblem = SWFlowProblem.SWFlowProblem(sw_model=0,
                                              cfl=0.33,
                                              outputStepping=outputStepping,
                                              structured=True,
                                              he=1.,
                                              nnx=1,
                                              nny=1,
                                              domain=domain,
                                              AdH_file=AdH_file,
                                              initialConditions=initialConditions,
                                              boundaryConditions=boundaryConditions,
                                              reflectingBCs=opts.reflectingBCs,
                                              bathymetry=bathymetry)
mySWFlowProblem.physical_parameters['LINEAR_FRICTION']=0
mySWFlowProblem.physical_parameters['mannings']=0.033
