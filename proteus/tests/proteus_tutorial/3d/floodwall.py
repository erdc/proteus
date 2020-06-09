"""
Multiphase Flow Test
"""
from __future__ import division
from past.utils import old_div
import numpy as np
from proteus import (Domain, Context, Gauges,
                     MeshTools as mt)
from proteus.Gauges import PointGauges, LineIntegralGauges, LineGauges
from proteus.Profiling import logEvent
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
import math
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus import WaveTools as wt
from proteus.mprans import SpatialTools as st

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ('ns_model',1,"ns_model={0,1} for {rans2p,rans3p}"),
    ("final_time",7.5,"Final time for simulation"),
    ("dt_output",0.1,"Time interval to output solution"),
    ("gauges", True, "Collect data for validation"),
    ("cfl",0.2,"Desired CFL restriction"),
    ("he",0.5,"Max mesh element diameter"),
    ("ARTIFICIAL_VISCOSITY",3,"artificial viscosity")
    ])

waterLevel=0.75
leeward_wl=0.25
pro_wl=0.5
g=np.array([0.,0.0,-9.81])
he=opts.he
# ****************** #
# ***** GAUGES ***** #
# ****************** #

# *************************** #
# ***** DOMAIN AND MESH ***** #
# *************************** #
domain = Domain.PiecewiseLinearComplexDomain()

boundaries=['gate','left','right','bottom','top','front','back']
boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])

boundaryOrientations = {'gate': np.array([-1., 0.,0.]),
                        'left': np.array([-1., 0.,0.]),
                        'right': np.array([+1., 0.,0.]),
                        'bottom': np.array([0., 0.,-1.]),
                        'top': np.array([0., 0.,+1.]),
                        'front': np.array([0.,+1.,0.]),
                        'back': np.array([0.,-1., 0.]),
                            }
vertices=[[0.0,0.0,0.0],#0
          [2.0,0.0,0.0],#1
          [2.0,0.0,1.0],#2
          [2.15,0.0,1.0],#3
          [2.75,0.0,0.0],#4
          [3.5,0.0,0.0],#5
          [3.5,0.0,1.5],#6
          [0.0,0.0,1.5],#7
          
          [0.0,2.,0.0],#8
          [2.0,2.,0.0],#9
          [2.0,2.,1.0],#10
          [2.15,2.,1.0],#11
          [2.75,2.,0.0],#12
          [3.5,2.,0.0],#13
          [3.5,2.,1.5],#14
          [0.0,2.,1.5],]#15

vertexFlags=[boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['top'],
             boundaryTags['top'],
             
             boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['bottom'],
             boundaryTags['top'],
             boundaryTags['top'],]

facets=[[[0,1,2,3,4,5,6,7]],
        [[8,9,10,11,12,13,14,15]],
        [[0,7,15,8]],
        [[13,14,6,5]],
        [[6,7,15,14]],
        [[0,1,9,8]],
        [[1,2,10,9]],
        [[2,3,11,10]],
        [[3,4,12,11]],
        [[4,5,13,12]]]



facetFlags=[boundaryTags['front'],
            boundaryTags['back'],
            boundaryTags['left'],
            boundaryTags['right'],
            boundaryTags['top'],
            boundaryTags['bottom'],
            boundaryTags['gate'],
            boundaryTags['gate'],
            boundaryTags['gate'],
            boundaryTags['bottom']]
        
regions=[[0.1,0.1,0.1]]
regionFlags=[1]


tank = st.CustomShape(domain,
                      vertices=vertices,
                      vertexFlags=vertexFlags,
                      facets=facets,
                      facetFlags=facetFlags,
                      regions = regions,
                      regionFlags = regionFlags,
                      boundaryTags=boundaryTags,
                      boundaryOrientations=boundaryOrientations)

tank.BC['top'].setAtmosphere()
#tank.BC['top'].setFreeSlip()
tank.BC['bottom'].setFreeSlip()
tank.BC['gate'].setFreeSlip()
tank.BC['front'].setFreeSlip()
tank.BC['back'].setFreeSlip()
#tank.BC['left'].setFreeSlip()
tank.BC['right'].setFreeSlip()


wave=wt.MonochromaticWaves(period=1.0,waveHeight=0.5,mwl=waterLevel,depth=waterLevel,g=g,waveDir=np.array([1.0,0.,0.]))
tank.BC['left'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=3.*he)

domain.MeshOptions.setParallelPartitioningType('node')
domain.boundaryTags = boundaryTags
he = opts.he
domain.MeshOptions.he = he
st.assembleDomain(domain)
domain.MeshOptions.triangleOptions="VApq1.25q12feena%e" % ((he**3)/6.0,)
domain.writePLY("mesh")
#domain.readPoly("SG")
domain.writePoly("mesh")
domain.writeAsymptote("mesh")

# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
class zero(object):
    def uOfXT(self,x,t):
        return 0.

def signedDistance(x):
    if x[0]>2.1:
        return x[2] - leeward_wl
    else:
        return x[2]-waterLevel

    # a = ((x[0]-19)**2+(x[1]-4.5)**2)**0.5
    # b = ((x[0]-19)**2+(x[1]+4.5)**2)**0.5

    # if x[0]>18.0:
    #     return x[2] - pro_wl
    # elif a<4 or b<4:
    #     return x[2] - pro_wl
    # else:
    #     return x[2] - waterLevel

epsFact_consrv_heaviside=3.0
class VF_IC:
    def uOfXT(self, x, t):
        return smoothedHeaviside(epsFact_consrv_heaviside*he,signedDistance(x))

class PHI_IC:
    def uOfXT(self, x, t):
        return signedDistance(x)


# ******************************* #
# ***** BOUNDARY CONDITIONS ***** #
# ******************************* #
    
# ADVECTIVE FLUX BOUNDARY CONDITIONS #
    
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
                                             nd=3,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             structured=False,
                                             he=he,
                                             nnx=None,
                                             nny=None,
                                             nnz=None,
                                             domain=domain,
                                             initialConditions=initialConditions,
                                             boundaryConditions=None,
                                             auxVariables=None,
                                             useSuperlu=False)
#myTpFlowProblem.physical_parameters['gravity'] = [0.0,0.0,-9.81]
myTpFlowProblem.Parameters.physical.gravity = [0., 0., -9.81]
