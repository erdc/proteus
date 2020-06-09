"""
dambreak 2-D
"""
from proteus import (Domain, Context)
from proteus.mprans.SpatialTools import Tank2D
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus.Gauges import PointGauges, LineIntegralGauges, LineGauges
import numpy as np


# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ("final_time",2.0,"Final time for simulation"),
    ("dt_output",0.01,"Time interval to output solution"),
    ("cfl",0.9,"Desired CFL restriction"),
    ("he",0.08,"he relative to Length of domain in x"),
    ("x_tank",3.22,"extent of domain in x"),
    ("y_tank",1.8,"extent of domain in y")
    ])

#Try adding these variables as Context options#
waterLine_y = 0.6 #the extent of the water column in +y from 0
waterLine_x = 1.2 #the extent of the water column in +x from 0


# *************************** #
# ***** DOMAIN AND MESH ***** #
# ****************** #******* #
tank_dim = (opts.x_tank,opts.y_tank) 

structured=False
if structured:
    domain = Domain.RectangularDomain(tank_dim)
else:
    domain = Domain.PlanarStraightLineGraphDomain()

# ----- TANK ----- #
tank = Tank2D(domain, tank_dim) 

# ----- BOUNDARY CONDITIONS ----- #
tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['x-'].setFreeSlip()


# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
class zero(object):
    def uOfXT(self,x,t):
        return 0.0

waterLine_y = 0.6
waterLine_x = 1.2

class PHI_IC:
    def uOfXT(self, x, t):
        phi_x = x[0] - waterLine_x
        phi_y = x[1] - waterLine_y
        if phi_x < 0.0:
            if phi_y < 0.0:
                return max(phi_x, phi_y)
            else:
                return phi_y
        else:
            if phi_y < 0.0:
                return phi_x
            else:
                return (phi_x ** 2 + phi_y ** 2)**0.5

class VF_IC:
    def __init__(self):
        self.phi = PHI_IC()
    def uOfXT(self, x, t):
        from proteus.ctransportCoefficients import smoothedHeaviside
        return smoothedHeaviside(1.5*opts.he,self.phi.uOfXT(x,0.0))

#########################
# ***** Numerics ****** #
#########################
#for structured
nny = int(tank_dim[1]/opts.he)+1
nnx = int(tank_dim[0]/opts.he)+1
#for unstructured
domain.MeshOptions.he = opts.he
st.assembleDomain(domain)
domain.MeshOptions.triangleOptions = "VApq30Dena%8.8f" % ((opts.he ** 2)/2.0,)

############################################
# ***** Create myTwoPhaseFlowProblem ***** #
############################################
outputStepping = TpFlow.OutputStepping(opts.final_time,dt_output=opts.dt_output,dt_init=0.0001)
initialConditions = {'pressure': zero(),
                     'pressure_increment': zero(),
                     'vel_u': zero(),
                     'vel_v': zero(),
                     'vof':  VF_IC(),
                     'ncls': PHI_IC(),
                     'rdls': PHI_IC()}

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=0,
                                             ls_model=0,
                                             nd=2,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             structured=structured,
                                             he=opts.he,
                                             nnx=nnx,
                                             nny=nny,
                                             domain=domain,
                                             initialConditions=initialConditions)
