"""
Rising Bubble Benchmarks
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
    ("test_case",1,"Rising bubble test cases"),
    ("final_time",3.0,"Final time for simulation"),
    ("dt_output",0.01,"Time interval to output solution"),
    ("cfl",0.33,"Desired CFL restriction"),
    ("refinement",3,"level of refinement"),
    ("structured",True,"use structured mesh")
    ])

# *************************** #
# ***** DOMAIN AND MESH ***** #
# *************************** #
tank_dim = (1.0,2.0)
# MESH
refinement = opts.refinement
structured = opts.structured
boundaries = ['bottom', 'right', 'top', 'left', 'front', 'back']
boundaryTags = dict([(key, i + 1) for (i, key) in enumerate(boundaries)])
if structured:
    nnx = 4 * refinement**2 + 2
    nny = 2 * nnx
    nnz = None
    domain = Domain.RectangularDomain(tank_dim)
    domain.boundaryTags = boundaryTags
    he = tank_dim[0]/(nnx - 1)
    triangleFlag=1
else:
    nnx = nny = None
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
    domain.boundaryTags = boundaryTags
    domain.writePoly("mesh")
    domain.writePLY("mesh")
    domain.writeAsymptote("mesh")
    he = old_div(tank_dim[0], float(4 * refinement - 1))
    domain.MeshOptions.he = he
    domain.MeshOptions.triangleOptions = "VApq30Dena%8.8f" % (old_div((he ** 2), 2.0),)

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
        # dist to center of bubble
        r = np.sqrt((x[0]-xB)**2 + (x[1]-yB)**2)
        # dist to surface of bubble
        dB = rB - r
        return dB

# ******************************* #
# ***** BOUNDARY CONDITIONS ***** #
# ******************************* #
# Non-slip BCs in top and bottom
# Slip BCs in left and right
# DIRICHLET BCs #
def vel_u_DBC(x,flag):
    if flag==boundaryTags['bottom'] or flag==boundaryTags['top']:
        return lambda x,t: 0.0

def vel_v_DBC(x,flag):
    if flag==boundaryTags['bottom'] or flag==boundaryTags['top']:
        return lambda x,t: 0.0

# ADVECTIVE FLUX BCs #
def vel_u_AFBC(x,flag):
    if flag==boundaryTags['left'] or flag==boundaryTags['right']:
        return lambda x,t: 0.

def vel_v_AFBC(x,flag):
    if flag==boundaryTags['left'] or flag==boundaryTags['right']:
        return lambda x,t: 0.

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
    'pressure_DBC': lambda x,t: None,
    'pressure_increment_DBC': lambda x,t: None,
    'vel_u_DBC': vel_u_DBC,
    'vel_v_DBC': vel_v_DBC,
    'clsvof_DBC': lambda x,t: None,
    # ADVECTIVE FLUX BCs #
    'pressure_AFBC': lambda x,flag: lambda x,t: 0,
    'pressure_increment_AFBC': lambda x,flag: lambda x,t: 0.,
    'vel_u_AFBC': vel_u_AFBC,
    'vel_v_AFBC': vel_v_AFBC,
    'clsvof_AFBC': lambda x,flag: lambda x,t: 0.,
    # DIFFUSIVE FLUX BCs #
    'pressure_increment_DFBC': lambda x,flag: lambda x,t: 0.,
    'vel_u_DFBC': lambda x,flag: lambda x,t: 0.,
    'vel_v_DFBC': lambda x,flag: lambda x,t: 0.,
    'clsvof_DFBC': lambda x,flag: lambda x,t: 0.}
myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=1,
                                             nd=2,
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
                                             useSuperlu=False)
physical_parameters = myTpFlowProblem.physical_parameters
physical_parameters['gravity'] = [0.0, -0.98, 0.0]
if opts.test_case==1:
    physical_parameters['densityA'] = 1000.0
    physical_parameters['viscosityA'] = 10.0/physical_parameters['densityA']
    physical_parameters['densityB'] = 100.0
    physical_parameters['viscosityB'] = 1.0/physical_parameters['densityB']
    physical_parameters['surf_tension_coeff'] = 24.5
    physical_parameters['gravity'] = [0.0, -0.98, 0.0]
else: #test_case=2
    physical_parameters['densityA'] = 1000.0
    physical_parameters['viscosityA'] = 10.0/physical_parameters['densityA']
    physical_parameters['densityB'] = 1.0
    physical_parameters['viscosityB'] = 0.1/physical_parameters['densityB']
    physical_parameters['surf_tension_coeff'] = 1.96
