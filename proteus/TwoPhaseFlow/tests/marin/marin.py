"""
Multiphase Flow Test
"""
from __future__ import division
from past.utils import old_div
import numpy as np
from proteus import (Domain, Context,
                     MeshTools as mt)
from proteus.Profiling import logEvent
from proteus.mprans.SpatialTools import Tank2D
from proteus.mprans.SpatialTools import Tank3D
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
import math

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ('ns_model',0,"ns_model = {rans2p,rans3p}"),
    ('ls_model',1,"ls_model = {ncls,clsvof}"),
    ("final_time",7.5,"Final time for simulation"),
    ("dt_output",0.01,"Time interval to output solution"),
    ("cfl",0.33,"Desired CFL restriction"),
    ("he",0.5,"Max mesh element diameter")
    ])

# ****************** #
# ***** GAUGES ***** #
# ****************** #
# None

# *************************** #
# ***** DOMAIN AND MESH ***** #
# ****************** #******* #
L      = [3.22,1.0,1.0]
box_L  = [0.161,0.403,0.161]
box_xy = [2.3955,0.2985]
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
regions=[[0.5*L[0],0.5*L[1],0.5*L[2]]]
regionFlags=[0]
domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                             vertexFlags=vertexFlags,
                                             facets=facets,
                                             facetFlags=facetFlags,
                                             regions = regions,
                                             regionFlags = regionFlags,
                                             holes=holes)
domain.MeshOptions.setParallelPartitioningType('node')
domain.boundaryTags = boundaryTags
domain.writePoly("mesh")
domain.writePLY("mesh")
domain.writeAsymptote("mesh")
domain.MeshOptions.triangleOptions="VApq1.25q12feena%e" % ((he**3)/6.0,)

# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
class zero(object):
    def uOfXT(self,x,t):
        return 0.
    
class clsvof_init_cond(object):
    def uOfXT(self,x,t):
        waterLine_x = 1.22
        waterLine_z = 0.55
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
                return math.sqrt(phi_x**2 + phi_z**2)

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
    
def pressure_increment_DBC(x,flag):
    if flag == boundaryTags['top'] and openTop:
        return lambda x,t: 0.0    

def pressure_DBC(x,flag):
    if flag == boundaryTags['top'] and openTop:
        return lambda x,t: 0.0

def clsvof_DBC(x,flag):
    if openTop and flag == boundaryTags['top']:
        return lambda x,t: 1.0
    
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

def pressure_increment_AFBC(x,flag):
    if not (flag == boundaryTags['top'] and openTop):
        return lambda x,t: 0.0

def pressure_AFBC(x,flag):
    if not(flag == boundaryTags['top'] and openTop):
        return lambda x,t: 0.0
    
def clsvof_AFBC(x,flag):
    if openTop and flag == boundaryTags['top']:
        return None
    else:
        return lambda x,t: 0.0

# DIFFUSIVE FLUX BCs #
def pressure_increment_DFBC(x,flag):
    if not (flag == boundaryTags['top'] and openTop):
        return lambda x,t: 0.0
    
############################################
# ***** Create myTwoPhaseFlowProblem ***** #
############################################
outputStepping = TpFlow.OutputStepping(opts.final_time,dt_output=opts.dt_output)
initialConditions = {'pressure': zero(),
                     'pressure_increment': zero(),
                     'vel_u': zero(),
                     'vel_v': zero(),
                     'vel_w': zero(),
                     'clsvof': clsvof_init_cond()}
boundaryConditions = {
    # DIRICHLET BCs #
    'pressure_DBC': pressure_DBC,
    'pressure_increment_DBC':  pressure_increment_DBC,
    'vel_u_DBC': vel_u_DBC,
    'vel_v_DBC': vel_v_DBC,
    'vel_w_DBC': vel_w_DBC,
    'clsvof_DBC': clsvof_DBC,
    # ADVECTIVE FLUX BCs #
    'pressure_AFBC': pressure_AFBC,
    'pressure_increment_AFBC': pressure_increment_AFBC,
    'vel_u_AFBC': vel_u_AFBC,
    'vel_v_AFBC': vel_v_AFBC,
    'vel_w_AFBC': vel_w_AFBC,
    'clsvof_AFBC': clsvof_AFBC,
    # DIFFUSIVE FLUX BCs #
    'pressure_increment_DFBC': pressure_increment_DFBC,
    'vel_u_DFBC': lambda x, flag: lambda x,t: 0.,
    'vel_v_DFBC': lambda x, flag: lambda x,t: 0.,
    'vel_w_DFBC': lambda x, flag: lambda x,t: 0.,
    'clsvof_DFBC': lambda x, flag: None}
myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=opts.ns_model,
                                             ls_model=opts.ls_model,
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
                                             boundaryConditions=boundaryConditions)
myTpFlowProblem.Parameters.physical['gravity'] = [0.0,0.0,-9.8]

params = myTpFlowProblem.Parameters

# MESH PARAMETERS
params.mesh.genMesh = opts.genMesh
params.mesh.he = he

# if opts.ns_model == 0:
#     myTpFlowProblem.Parameters.Models.rans2p['index'] = 0
#     myTpFlowProblem.Parameters.Models.clsvof['index'] = 1
# elif opts.ns_model == 1:
#     myTpFlowProblem.Parameters.Models.clsvof['index'] = 0
#     myTpFlowProblem.Parameters.Models.rans3p['index'] = 1
#     myTpFlowProblem.Parameters.Models.pressureIncrement['index'] = 2
#     myTpFlowProblem.Parameters.Models.pressure['index'] = 3
#     myTpFlowProblem.Parameters.Models.pressureInitial['index'] = 4
