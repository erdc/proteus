"""
Multiphase Flow Test
"""
from __future__ import division
from past.utils import old_div
import numpy as np
from proteus import (Domain, Context,
                     MeshTools as mt)
from proteus.Profiling import logEvent
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ('nd',2,"Num of dimensions"),
    ('ns_model',0,"ns_model = {rans2p,rans3p}"),
    ("final_time",3.0,"Final time for simulation"),
    ("dt_output",0.01,"Time interval to output solution"),
    ("cfl",0.33,"Desired CFL restriction"),
    ("refinement",3,"level of refinement")
    ],mutable=True)

# ****************** #
# ***** GAUGES ***** #
# ****************** #
# None

# *************************** #
# ***** DOMAIN AND MESH ***** #
# ****************** #******* #
tank_dim = (1.0,1.0) if opts.nd == 2 else (1.0,1.0,1.0)
# MESH
refinement = opts.refinement
structured=True
boundaries = ['bottom', 'right', 'top', 'left', 'front', 'back']
boundaryTags = dict([(key, i + 1) for (i, key) in enumerate(boundaries)])
if structured:
    nnx = 4 * refinement**2 + 1
    nny = nnx
    nnz = None if opts.nd == 2 else nnx
    domain = Domain.RectangularDomain(tank_dim)
    domain.boundaryTags = boundaryTags
    he = tank_dim[0]/(nnx - 1)
else:
    if opts.nd==2:
        vertices = [[0.0, 0.0],  #0
                    [tank_dim[0], 0.0],  #1
                    [tank_dim[0], tank_dim[1]],  #2
                    [0.0, tank_dim[1]]]  #3
        vertexFlags = [boundaryTags['bottom'],
                       boundaryTags['bottom'],
                       boundaryTags['top'],
                       boundaryTags['top']]
        segments = [[0, 1],
                    [1, 2],
                    [2, 3],
                    [3, 0]]
        segmentFlags = [boundaryTags['bottom'],
                        boundaryTags['right'],
                        boundaryTags['top'],
                        boundaryTags['left']]
        regions = [[tank_dim[0]/2., tank_dim[1]/2.]]
        regionFlags = [1]
        domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                      vertexFlags=vertexFlags,
                                                      segments=segments,
                                                      segmentFlags=segmentFlags,
                                                      regions=regions,
                                                      regionFlags=regionFlags)
    else: #nd==3
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
        regions=[[0.5*L[0],0.5*L[1],0.5*L[2]]]
        regionFlags=[1]
        domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                     vertexFlags=vertexFlags,
                                                     facets=facets,
                                                     facetFlags=facetFlags,
                                                     regions=regions,
                                                     regionFlags=regionFlags)
    domain.boundaryTags = boundaryTags
    domain.writePoly("mesh")
    domain.writePLY("mesh")
    domain.writeAsymptote("mesh")
    he = old_div(tank_dim[0], float(4 * refinement - 1))
    domain.MeshOptions.he = he
    if opts.nd==2:
        triangleOptions = "VApq30Dena%8.8f" % (old_div((he ** 2), 2.0),)
    else:
        triangleOptions="VApq1.4q12feena%21.16e" % (old_div((he**3),6.0),)

# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
class zero(object):
    def uOfXT(self,x,t):
        return 0.
class clsvof_init_cond(object):
    def uOfXT(self,x,t):
        xB = 0.5
        yB = 0.5
        rB = 0.25
        zB = 0.5
        # dist to center of bubble
        if opts.nd==2:
            r = np.sqrt((x[0]-xB)**2 + (x[1]-yB)**2)
        else:
            r = np.sqrt((x[0]-xB)**2 + (x[1]-yB)**2 + (x[2]-zB)**2)
        # dist to surface of bubble
        dB = -(rB - r)
        return dB

# ******************************* #
# ***** BOUNDARY CONDITIONS ***** #
# ******************************* #
#################
# NAVIER STOKES #
#################
# DIRICHLET BCs #
def pressure_DBC(x,flag):
    if (flag==boundaryTags['top']):
        return lambda x,t: 0.
def pressure_increment_DBC(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0
vel_u_DBC = lambda x,flag: None
vel_v_DBC = lambda x,flag: None
vel_w_DBC = lambda x,flag: None

# ADVECTIVE FLUX #
def pressure_AFBC(x,flag):
    if not (flag==boundaryTags['top']):
        return lambda x,t: 0.
def pressure_increment_AFBC(x,flag):
    if not (flag == boundaryTags['top']):
        return lambda x,t: 0.0
def vel_u_AFBC(x,flag):
    if not (flag==boundaryTags['top']):
        return lambda x,t: 0.
def vel_v_AFBC(x,flag):
    if not (flag==boundaryTags['top']):
        return lambda x,t: 0.
def vel_w_AFBC(x,flag):
    if not (flag==boundaryTags['top']):
        return lambda x,t: 0.

# DIFFUSIVE FLUX #
def pressure_increment_DFBC(x,flag):
    if not (flag == boundaryTags['top']):
        return lambda x,t: 0.0
vel_u_DFBC = lambda x,flag: lambda x,t: 0.0
vel_v_DFBC = lambda x,flag: lambda x,t: 0.0
vel_w_DFBC = lambda x,flag: lambda x,t: 0.0

#############
# LEVEL SET #
#############
def clsvof_DBC(x,flag):
    if (flag==boundaryTags['top']): #Just let air in
        return lambda x,t: 1.
def clsvof_AFBC(x,flag):
    if not (flag==boundaryTags['top']):
        return lambda x,t: 0.
clsvof_DFBC = lambda x,flag: lambda x,t: 0.0

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
    'pressure_increment_DBC': pressure_increment_DBC,
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
    'vel_u_DFBC': vel_u_DFBC,
    'vel_v_DFBC': vel_v_DFBC,
    'vel_w_DFBC': vel_w_DFBC,
    'clsvof_DFBC': clsvof_DFBC}
myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=opts.ns_model,
                                             nd=opts.nd,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             structured=structured,
                                             he=he,
                                             nnx=nnx,
                                             nny=nny,
                                             nnz=None,
                                             domain=domain,
                                             initialConditions=initialConditions,
                                             boundaryConditions=boundaryConditions,
                                             useSuperlu=True)
myTpFlowProblem.useBoundaryConditionsModule = False
myTpFlowProblem.Parameters.Models.rans2p.epsFact_viscosity = 3.
myTpFlowProblem.Parameters.Models.rans2p.epsFact_density = 3.
myTpFlowProblem.Parameters.Models.rans2p.ns_shockCapturingFactor = 0.25
myTpFlowProblem.Parameters.Models.rans2p.timeDiscretization = 'vbdf'

myTpFlowProblem.Parameters.physical.gravity = np.array([0., -9.8, 0.])

myTpFlowProblem.outputStepping.systemStepExact = True
