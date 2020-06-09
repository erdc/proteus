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
from subprocess import check_call
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus.TwoPhaseFlow.utils.Parameters import ParametersPhysical as PP
from proteus.ctransportCoefficients import smoothedHeaviside
import math

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ("final_time",7.5,"Final time for simulation"),
    ("dt_output",0.1,"Time interval to output solution"),
    ("gauges", True, "Collect data for validation"),
    ("cfl",0.2,"Desired CFL restriction"),
    ("he",0.5,"Max mesh element diameter"),
    ("use_gmsh",False,"Use gmsh to generate mesh")
    ])

#assert opts.ns_model==1, "use ns_model=1 (rans3pf) for this"

# ****************** #
# ***** GAUGES ***** #
# ****************** #
if opts.gauges:
    pressure_gauges = PointGauges(gauges=((('p',),
                                           ((2.389,0.526,0.025), #P1
                                            (2.389,0.526,0.099), #P3
                                            (2.414,0.474,0.165), #P5
                                            (2.487,0.474,0.165))),), #P7
                                  fileName="pressure.csv")
    point_height_gauges = PointGauges(gauges=((('phi',),
                                               ((2.389,0.526,0.025), #P1
                                                (2.389,0.526,0.099), #P3
                                                (2.414,0.474,0.165), #P5
                                                (2.487,0.474,0.165))),), #P7
                                      fileName="point_clsvof.csv")
    height_gauges = LineGauges(gauges=((("phi",),
                                        (((2.724, 0.5, 0.0),
                                          (2.724, 0.5, 1.0)),
                                         ((2.228, 0.5, 0.0),
                                          (2.228, 0.5, 1.0)),
                                         ((1.732, 0.5, 0.0),
                                          (1.732, 0.5, 1.0)),
                                         ((0.582, 0.5, 0.0),
                                          (0.582, 0.5, 1.0)))),),
                               fileName="height.csv")

# *************************** #
# ***** DOMAIN AND MESH ***** #
# ****************** #******* #
L      = [3.22,1.0,1.0]
box_L  = [0.16,0.4,0.16]
box_xy = [2.39,0.3]
he = opts.he
boundaries=['left','right','bottom','top','front','back','box_left','box_right','box_top','box_front','box_back',]
boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
bt = boundaryTags
holes = [[0.5*box_L[0]+box_xy[0],0.5*box_L[1]+box_xy[1],0.5*box_L[2]]]
vertices=[[0.0,0.0,0.0],#0
          [L[0],0.0,0.0],#1
          [L[0],L[1],0.0],#2
          [0.0,L[1],0.0],#3
          [0.0,0.0,L[2]],#4
          [L[0],0.0,L[2]],#5
          [L[0],L[1],L[2]],#6
          [0.0,L[1],L[2]],#7
          [box_xy[0],box_xy[1],0.0],#8
          [box_xy[0]+box_L[0],box_xy[1],0.0],#9
          [box_xy[0]+box_L[0],box_xy[1]+box_L[1],0.0],#10
          [box_xy[0],box_xy[1]+box_L[1],0.0],#11
          [box_xy[0],box_xy[1],box_L[2]],#12
          [box_xy[0]+box_L[0],box_xy[1],box_L[2]],#13
          [box_xy[0]+box_L[0],box_xy[1]+box_L[1],box_L[2]],#14
          [box_xy[0],box_xy[1]+box_L[1],box_L[2]]]#15
vertexFlags=[boundaryTags['left'],
             boundaryTags['right'],
             boundaryTags['right'],
             boundaryTags['left'],
             boundaryTags['left'],
             boundaryTags['right'],
             boundaryTags['right'],
             boundaryTags['left'],
             boundaryTags['box_left'],
             boundaryTags['box_left'],
             boundaryTags['box_left'],
             boundaryTags['box_left'],
             boundaryTags['box_left'],
             boundaryTags['box_left'],
             boundaryTags['box_left'],
             boundaryTags['box_left']]
facets=[[[0,1,2,3],[8,9,10,11]],
        [[0,1,5,4]],
        [[1,2,6,5]],
        [[2,3,7,6]],
        [[3,0,4,7]],
        [[4,5,6,7]],
        [[8,9,13,12]],
        [[9,10,14,13]],
        [[10,11,15,14]],
        [[11,8,12,15]],
        [[12,13,14,15]]]
facetFlags=[boundaryTags['bottom'],
            boundaryTags['front'],
            boundaryTags['right'],
            boundaryTags['back'],
            boundaryTags['left'],
            boundaryTags['top'],
            boundaryTags['box_front'],
            boundaryTags['box_right'],
            boundaryTags['box_back'],
            boundaryTags['box_left'],
            boundaryTags['box_top']]
volumes=[[[fi for fi in range(len(facetFlags))]]]
tank_only=False
if tank_only:
    vertices=[[0.0,0.0,0.0],#0
              [L[0],0.0,0.0],#1
              [L[0],L[1],0.0],#2
              [0.0,L[1],0.0],#3
              [0.0,0.0,L[2]],#4
              [L[0],0.0,L[2]],#5
              [L[0],L[1],L[2]],#6
              [0.0,L[1],L[2]]]#7
    vertexFlags=[boundaryTags['left'],
                 boundaryTags['right'],
                 boundaryTags['right'],
                 boundaryTags['left'],
                 boundaryTags['left'],
                 boundaryTags['right'],
                 boundaryTags['right'],
                 boundaryTags['left']]
    facets=[[[0,1,2,3]],
            [[0,1,5,4]],
            [[1,2,6,5]],
            [[2,3,7,6]],
            [[3,0,4,7]],
            [[4,5,6,7]]]
    facetFlags=[boundaryTags['bottom'],
                boundaryTags['front'],
                boundaryTags['right'],
                boundaryTags['back'],
                boundaryTags['left'],
                boundaryTags['top']]
    volumes=[[[0,1,2,3,4,5]]]
regions=[[0.5*L[0],0.5*L[1],0.5*L[2]]]
regionFlags=[0]
domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                             vertexFlags=vertexFlags,
                                             facets=facets,
                                             facetFlags=facetFlags,
                                             regions = regions,
                                             regionFlags = regionFlags)
#                                             holes=holes)
domain.holes_ind=[]
domain.volumes=volumes
domain.MeshOptions.setParallelPartitioningType('node')
#domain.MeshOptions.use_gmsh=True
domain.boundaryTags = boundaryTags
domain.writePoly("mesh")
from proteus import Comm
if opts.use_gmsh:
    comm = Comm.get()
    if comm.isMaster():
        domain.writeGeo("mesh",he_max=he)
        domain.MeshOptions.triangleOptions="VApq1.25q12feena%e" % ((he**3)/6.0,)
        gmsh_cmd = "gmsh {0:s} -v 10 -3 -o {1:s} -format msh2".format(domain.geofile+".geo", domain.geofile+".msh")
        check_call(gmsh_cmd, shell=True)
        mt.msh2simplex("mesh",nd=3)
domain.MeshOptions.genMesh=False
# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
class zero(object):
    def uOfXT(self,x,t):
        return 0.

# Initial condition
waterLine_x = 1.20
waterLine_z = 0.55

def signedDistance(x):
    from math import sqrt
    phi_x = x[0]-waterLine_x
    phi_z = x[2]-waterLine_z 
    if phi_x < 0.0:
        if phi_z < 0.0:
            return max(phi_x,phi_z)
        else:
            return phi_z
    else:
        if phi_z < 0.0:
            return phi_x
        else:
            return sqrt(phi_x**2 + phi_z**2)

pp = PP()
rho_1 = pp.densityA
rho_0 = pp.densityB
g = [0., 0., -9.81]
class P:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        if signedDistance(x) < 0:
            return -(L[2] - self.waterLevel)*rho_1*g[2] - (self.waterLevel - x[2])*rho_0*g[2]
        else:
            return -(L[2] - x[2])*rho_1*g[2]

class V:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Phi:       
    def uOfXT(self,x,t):
        return signedDistance(x)

class VOF:
    def uOfXT(self,x,t):
        return smoothedHeaviside(1.5*he,signedDistance(x))



# ******************************* #
# ***** BOUNDARY CONDITIONS ***** #
# ******************************* #
non_slip_BCs=True
openTop=True
# DIRICHLET BOUNDARY CONDITIONS #
def vel_u_DBC(x,flag):
    if non_slip_BCs and (flag == boundaryTags['box_left'] or
                         flag == boundaryTags['box_right'] or
                         flag == boundaryTags['box_top'] or
                         flag == boundaryTags['box_front'] or
                         flag == boundaryTags['box_back']):
        return lambda  x,t: 0.0

def vel_v_DBC(x,flag):
    if non_slip_BCs and (flag == boundaryTags['box_left'] or
                         flag == boundaryTags['box_right'] or
                         flag == boundaryTags['box_top'] or
                         flag == boundaryTags['box_front'] or
                         flag == boundaryTags['box_back']):
        return lambda  x,t: 0.0

def vel_w_DBC(x,flag):
    if non_slip_BCs and (flag == boundaryTags['box_left'] or
                         flag == boundaryTags['box_right'] or
                         flag == boundaryTags['box_top'] or
                         flag == boundaryTags['box_front'] or
                         flag == boundaryTags['box_back']):
        return lambda  x,t: 0.0
    
def pressure_DBC(x,flag):
    if flag == boundaryTags['top'] and openTop:
        return lambda x,t: 0.0

def vof_DBC(x,flag):
    if openTop and flag == boundaryTags['top']:
        return lambda x,t: 1.0

def phi_DBC(x,flag):
    return None
    
# ADVECTIVE FLUX BOUNDARY CONDITIONS #
def vel_u_AFBC(x,flag):
    if non_slip_BCs and (flag == boundaryTags['box_left'] or
                         flag == boundaryTags['box_right'] or
                         flag == boundaryTags['box_top'] or
                         flag == boundaryTags['box_front'] or
                         flag == boundaryTags['box_back']):
        return None
    elif openTop and flag == boundaryTags['top']:
        return None
    else: #slip everywhere but the box
        return lambda x,t: 0.0

def vel_v_AFBC(x,flag):
    if non_slip_BCs and (flag == boundaryTags['box_left'] or
                         flag == boundaryTags['box_right'] or
                         flag == boundaryTags['box_top'] or
                         flag == boundaryTags['box_front'] or
                         flag == boundaryTags['box_back']):
        return None
    elif openTop and flag == boundaryTags['top']:
        return None
    else: #slip everywhere but the box
        return lambda x,t: 0.0

def vel_w_AFBC(x,flag):
    if non_slip_BCs and (flag == boundaryTags['box_left'] or
                         flag == boundaryTags['box_right'] or
                         flag == boundaryTags['box_top'] or
                         flag == boundaryTags['box_front'] or
                         flag == boundaryTags['box_back']):
        return None
    elif openTop and flag == boundaryTags['top']:
        return None
    else: #slip everywhere but the box
        return lambda x,t: 0.0

def pressure_AFBC(x,flag):
    if not(flag == boundaryTags['top'] and openTop):
        return lambda x,t: 0.0
    
def vof_AFBC(x,flag):
    if openTop and flag == boundaryTags['top']:
        return None
    else:
        return lambda x,t: 0.0

    
############################################
# ***** Create myTwoPhaseFlowProblem ***** #
############################################
outputStepping = TpFlow.OutputStepping(opts.final_time,dt_output=opts.dt_output)
initialConditions = {'pressure': P(waterLine_z),
                     'vel_u': zero(),
                     'vel_v': zero(),
                     'vel_w': zero(),
                     'vof': VOF(),
                     'rdls': Phi(),
                     'ncls': Phi()}

boundaryConditions = {
    # DIRICHLET BCs #
    'pressure_DBC': pressure_DBC,
    'vel_u_DBC': vel_u_DBC,
    'vel_v_DBC': vel_v_DBC,
    'vel_w_DBC': vel_w_DBC,
    'vof_DBC':  vof_DBC,
    'ncls_DBC': phi_DBC,
    # ADVECTIVE FLUX BCs #
    'pressure_AFBC': pressure_AFBC,
    'vel_u_AFBC': vel_u_AFBC,
    'vel_v_AFBC': vel_v_AFBC,
    'vel_w_AFBC': vel_w_AFBC,
    'vof_AFBC': vof_AFBC,
    # DIFFUSIVE FLUX BCs #
    'vel_u_DFBC': lambda x, flag: lambda x,t: 0.,
    'vel_v_DFBC': lambda x, flag: lambda x,t: 0.,
    'vel_w_DFBC': lambda x, flag: lambda x,t: 0.,
    'vof_DFBC': lambda x, flag: None}

auxVariables={'vof': [point_height_gauges, height_gauges],
              'pressure': [pressure_gauges]}

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=0,
                                             ls_model=0,
                                             nd=3,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             he=he,
                                             domain=domain,
                                             initialConditions=initialConditions,
                                             boundaryConditions=boundaryConditions,
                                             auxVariables=auxVariables,
                                             useSuperlu=True)
myTpFlowProblem.Parameters.physical.gravity = g
