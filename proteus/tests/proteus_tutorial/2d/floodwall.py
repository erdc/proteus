"""
Multiphase Flow Test
"""
import numpy as np
from proteus import (Domain, Context)
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
import math
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus import WaveTools as wt
from proteus.mprans import SpatialTools as st

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ('ns_model',1,"ns_model={0,1} for {rans2p,rans3p}"),
    ("final_time",10.0,"Final time for simulation"),
    ("dt_output",0.1,"Time interval to output solution"),
    ("gauges", True, "Collect data for validation"),
    ("cfl",0.9,"Desired CFL restriction"),
    ("he",0.08,"Max mesh element diameter")
    ])

waterLevel=0.9
pro_wl=0.5
g=np.array([0.,-9.81,0.0])
he=opts.he
# ****************** #
# ***** GAUGES ***** #
# ****************** #

# *************************** #
# ***** DOMAIN AND MESH ***** #
# *************************** #
domain = Domain.PlanarStraightLineGraphDomain()

boundaries=['gate','left','right','bottom','top']
boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])

boundaryOrientations = {'gate': np.array([0, 1.,0.]),
                        'left': np.array([-1., 0.,0.]),
                        'right': np.array([+1., 0.,0.]),
                        'top': np.array([0.,+1.,0.]),
                        'bottom': np.array([0.,-1., 0.]),
}
vertices=[[0.0,0.0],#0
          [2.0,0.0],#1
          [2.1,1.0],#2
          [2.15,1.0],#3
#          [2.75,0.0],#4
#          [3.5,0.0],#5
          [2.75 - 0.75*(2.75-2.15)/(1.0-0.0),0.75],#4
#          [3.5,0.75],#5
#          [3.5,1.5],#6
          [3.,0.75],#5
          [3.,1.5],#6
          [0.0,1.5],#7
]

vertexFlags=[boundaryTags['bottom'],
             boundaryTags['gate'],
             boundaryTags['gate'],
             boundaryTags['gate'],
             boundaryTags['gate'],
             boundaryTags['bottom'],
             boundaryTags['top'],
             boundaryTags['top'],
             ]

segments=[[0,1],
          [1,2],
          [2,3],
          [3,4],
          [4,5],
          [5,6],
          [6,7],
          [7,0]]



segmentFlags=[boundaryTags['bottom'],
              boundaryTags['gate'],
              boundaryTags['gate'],
              boundaryTags['gate'],
              boundaryTags['bottom'],
              boundaryTags['right'],
              boundaryTags['top'],
              boundaryTags['left']]

regions=[[0.1,0.1]]
regionFlags=[1]

tank = st.CustomShape(domain,
                      vertices=vertices,
                      vertexFlags=vertexFlags,
                      segments=segments,
                      segmentFlags=segmentFlags,
                      regions = regions,
                      regionFlags = regionFlags,
                      boundaryTags=boundaryTags,
                      boundaryOrientations=boundaryOrientations)

# ******************************* #
# ***** BOUNDARY CONDITIONS ***** #
# ******************************* #

wave=wt.MonochromaticWaves(period=2.5,
                           waveHeight=waterLevel/3.0,
                           mwl=waterLevel,
                           depth=waterLevel,
                           g=g,
                           waveDir=np.array([1.,0.,0.]))
tank.BC['left'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=1.5*he)
tank.BC['bottom'].setFreeSlip()
tank.BC['gate'].setFreeSlip()
tank.BC['right'].setFreeSlip()
tank.BC['top'].setAtmosphere()

domain.MeshOptions.setParallelPartitioningType('node')
domain.boundaryTags = boundaryTags
domain.MeshOptions.he = opts.he
st.assembleDomain(domain)
triangleOptions = "VApq30Dena%8.8f" % ((opts.he**2)/2.0,)   
domain.writePoly("mesh")

# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
class zero(object):
    def uOfXT(self,x,t):
        return 0.

def signedDistance(x):
    phi_x = x[0] - 2.15
    phi_y = x[1] - waterLevel
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
    def uOfXT(self, x, t):
        return smoothedHeaviside(1.5*opts.he,signedDistance(x))

class PHI_IC:
    def uOfXT(self, x, t):
        return signedDistance(x)
    
############################################
# ***** Create myTwoPhaseFlowProblem ***** #
############################################
outputStepping = TpFlow.OutputStepping(opts.final_time,dt_output=opts.dt_output)
initialConditions = {'pressure': zero(),
                     'pressure_increment': zero(),
                     'vel_u': zero(),
                     'vel_v': zero(),
                     'vel_w': zero(),
                     'vof': VF_IC(),
                     'ncls': PHI_IC(),
                     'rdls': PHI_IC(),}

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=0,
                                             ls_model=0,
                                             nd=2,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             he=he,
                                             domain=domain,
                                             initialConditions=initialConditions)
