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

waterLevel=2
pro_wl=0.5
g=np.array([0.,0.,-9.81])
he=opts.he
# ****************** #
# ***** GAUGES ***** #
# ****************** #

# *************************** #
# ***** DOMAIN AND MESH ***** #
# *************************** #
from proteus.mprans import SpatialTools as st
domain = Domain.PiecewiseLinearComplexDomain()
domain2 = Domain.PiecewiseLinearComplexDomain()

SG=st.ShapeSTL(domain2,'SG_full_2.stl')


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


facetFlags=SG.facetFlags.tolist()
vertexFlags=SG.vertexFlags.tolist()
facets=SG.facets.tolist()
vertices=SG.vertices.tolist()

for i in range(len(SG.vertexFlags.tolist())):
    if SG.vertices.tolist()[i][0]==min(SG.vertices[:,0]):
        vertexFlags[i]=2
    elif SG.vertices.tolist()[i][0]==max(SG.vertices[:,0]):
        vertexFlags[i]=3
    elif SG.vertices.tolist()[i][1]==min(SG.vertices[:,1]):
        vertexFlags[i]=6
    elif SG.vertices.tolist()[i][1]==max(SG.vertices[:,1]):
        vertexFlags[i]=7
    elif SG.vertices.tolist()[i][2]==min(SG.vertices[:,2]):
        vertexFlags[i]=4
    elif SG.vertices.tolist()[i][2]==max(SG.vertices[:,2]):
        vertexFlags[i]=5
    else:
        vertexFlags[i]=1
for i in range(len(facetFlags)):
    if (vertices[facets[i][0][0]][0]==min(SG.vertices[:,0])
        and vertices[facets[i][0][1]][0]==min(SG.vertices[:,0])
        and vertices[facets[i][0][2]][0]==min(SG.vertices[:,0])):
        facetFlags[i]=2
    elif (vertices[facets[i][0][0]][0]==max(SG.vertices[:,0])
          and vertices[facets[i][0][1]][0]==max(SG.vertices[:,0])
          and vertices[facets[i][0][2]][0]==max(SG.vertices[:,0])):
        facetFlags[i]=3
    elif (vertices[facets[i][0][0]][1]==min(SG.vertices[:,1])
          and vertices[facets[i][0][1]][1]==min(SG.vertices[:,1])
          and vertices[facets[i][0][2]][1]==min(SG.vertices[:,1])):
        facetFlags[i]=6
    elif (vertices[facets[i][0][0]][1]==max(SG.vertices[:,1])
          and vertices[facets[i][0][1]][1]==max(SG.vertices[:,1])
          and vertices[facets[i][0][2]][1]==max(SG.vertices[:,1])):
        facetFlags[i]=7
    elif(vertices[facets[i][0][0]][2]==min(SG.vertices[:,2])
         and vertices[facets[i][0][1]][2]==min(SG.vertices[:,2])
         and vertices[facets[i][0][2]][2]==min(SG.vertices[:,2])):
        facetFlags[i]=4
    elif (vertices[facets[i][0][0]][2]==max(SG.vertices[:,2])
          and vertices[facets[i][0][1]][2]==max(SG.vertices[:,2])
          and vertices[facets[i][0][2]][2]==max(SG.vertices[:,2])):
        facetFlags[i]=5
    else:
        facetFlags[i]=1
regions=[[0.5,0.5,0.5]]
regionFlags=[1]
holes=[[18.0,0.0,1.5]]
L=[max(SG.vertices[:,0]),max(SG.vertices[:,1]),max(SG.vertices[:,2])]
tank = st.CustomShape(domain,
                      vertices=vertices,
                      vertexFlags=vertexFlags,
                      facets=facets,
                      facetFlags=facetFlags,
                      regions = regions,
                      regionFlags = regionFlags,
                      holes=holes,
                      boundaryTags=boundaryTags,
                      boundaryOrientations=boundaryOrientations)

tank.BC['top'].setAtmosphere()
tank.BC['bottom'].setFreeSlip()
tank.BC['gate'].setFreeSlip()
tank.BC['front'].setFreeSlip()
tank.BC['back'].setFreeSlip()
tank.BC['right'].setFreeSlip()


wave=wt.MonochromaticWaves(period=2.0,waveHeight=1.0,mwl=waterLevel,depth=waterLevel,g=g,waveDir=np.array([1.0,0.,0.]))
tank.BC['left'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=3.*he)

domain.MeshOptions.setParallelPartitioningType('node')
domain.boundaryTags = boundaryTags
he = opts.he
domain.MeshOptions.he = he
st.assembleDomain(domain)
domain.MeshOptions.triangleOptions="VApq1.25q12feena%e" % ((he**3)/6.0,)
domain.writePLY("mesh")
domain.writePoly("mesh")
domain.writeAsymptote("mesh")

# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
class zero(object):
    def uOfXT(self,x,t):
        return 0.

def signedDistance(x):
    a = ((x[0]-19.5)**2+(x[1]-4.5)**2)**0.5
    b = ((x[0]-19.5)**2+(x[1]+4.5)**2)**0.5

    if x[0]>18.0:
        return x[2] - pro_wl
    elif a<4 or b<4:
        return x[2] - pro_wl
    else:
        return x[2] - waterLevel

epsFact_consrv_heaviside=3.0
class VF_IC:
    def uOfXT(self, x, t):
        return smoothedHeaviside(epsFact_consrv_heaviside*he,signedDistance(x))

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
                                             nd=3,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             he=he,
                                             domain=domain,
                                             initialConditions=initialConditions,
                                             useSuperlu=False)
myTpFlowProblem.Parameters.physical.gravity = [0., 0., -9.81]
