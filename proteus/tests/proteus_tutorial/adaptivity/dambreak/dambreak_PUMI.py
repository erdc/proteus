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
he = opts.he

#initialize PUMI domain
domain = Domain.PUMIDomain(dim=2) #initialize the domain

#read the geometry and mesh
#domain.AdaptManager.PUMIAdapter.loadModelAndMesh(b"Reconstructed.dmg", b"Reconstructed.smb")
domain.AdaptManager.PUMIAdapter.loadModelAndMesh(b"Reconstructed.dmg", b"4-Proc/.smb")

#mesh is partitioned by elements and not by nodes
from proteus import MeshTools
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.element
domain.MeshOptions.setParallelPartitioningType('element')

#input for adaptation describing the staggered solve, needed for proper post-adapt progression
domain.AdaptManager.modelDict = {'flow':0,'phase':2,'correction':[3,4]}

#input to determine what strategy for adaptation is used
domain.AdaptManager.sizeInputs = [b'interface',b'error_vms']

#do you want to adapt?
domain.AdaptManager.adapt = 1

#max,min mesh edge length
domain.AdaptManager.hmax = he*2.0
domain.AdaptManager.hmin= he/2.0

#mesh edge length for near interface
domain.AdaptManager.hphi= he/2.0

#number of predictive timesteps for interface-based adapt
domain.AdaptManager.numAdaptSteps= 10

#target element error for error estimation-based adapt
domain.AdaptManager.targetError= 2.0

#maximum ratio between adjacent edge lengths
domain.AdaptManager.gradingFactor= 1.5

#number of iterations in adapt routine
domain.AdaptManager.numIterations= 5

#see before and after meshes and other info for adaptation
domain.AdaptManager.logging= 0

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

#class PHI_IC:
#    def uOfXT(self, x, t):
#        phi_x = x[0] - waterLine_x
#        phi_y = x[1] - waterLine_y
#        if phi_x < 0.0:
#            if phi_y < 0.0:
#                return max(phi_x, phi_y)
#            else:
#                return phi_y
#        else:
#            if phi_y < 0.0:
#                return phi_x
#            else:
#                return (phi_x ** 2 + phi_y ** 2)**0.5

class PHI_IC:
    def uOfXT(self, x, t):
        import math
        radius = 0.5*waterLine_y
        cx = waterLine_x
        cy = waterLine_y
        R = math.sqrt((x[0] - cx)**2 + (x[1] - cy)**2)
        return R-radius

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

myTpFlowProblem.Parameters.Models.rans2p.n.conservativeFlux= None
