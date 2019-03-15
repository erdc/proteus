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
    ('ns_model',1,"ns_model = {rans2p,rans3p}"),
    ("final_time",3.0,"Final time for simulation"),
    ("dt_output",0.01,"Time interval to output solution"),
    ("cfl",0.33,"Desired CFL restriction"),
    ("refinement",16,"level of refinement")
    ],mutable=True)

# ****************** #
# ***** GAUGES ***** #
# ****************** #
# None

# *************************** #
# ***** DOMAIN AND MESH ***** #
# ****************** #******* #
tank_dim = (1.0,1.0)
# MESH
refinement=opts.refinement #64 if structured=False
structured=False
boundaries = ['bottom', 'right', 'top', 'left', 'front', 'back']
boundaryTags = dict([(key, i + 1) for (i, key) in enumerate(boundaries)])
if structured:
    nnx = 4 * refinement**2 + 1
    nny = nnx
    nnz = None
    domain = Domain.RectangularDomain(tank_dim)
    domain.boundaryTags = boundaryTags
    he = tank_dim[0]/(nnx - 1)
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
    def uOfXT(self,X,t):
        flat_inlet = False
        x=X[0]
        y=X[1]
        if flat_inlet:
            if x-tank_dim[0]/2. >= -0.05 and x-tank_dim[0]/2. <= 0.05:
                if y==tank_dim[1]:
                    return 0.
                else:
                    return tank_dim[1] - y
            else:
                return min(np.sqrt( (x-(tank_dim[0]/2.-0.05))**2 + (y-tank_dim[1])**2),
                           np.sqrt( (x-(tank_dim[0]/2.+0.05))**2 + (y-tank_dim[1])**2) )
        else: #circular inlet
            rInflow = np.sqrt((x-0.5)**2 + (y-1.0)**2)
            dB = -(0.05-rInflow)
            return dB

# ******************************* #
# ***** BOUNDARY CONDITIONS ***** #
# ******************************* #
#################
# NAVIER STOKES #
#################
# DIRICHLET BCs #
def pressure_DBC(x,flag):
    if flag == boundaryTags['top'] and not (x[0]>=0.5-0.05 and x[0]<=0.5+0.05):
        return lambda x,t: 0.0

def pressure_increment_DBC(x,flag):
    if flag == boundaryTags['top'] and not(x[0]>=0.5-0.05 and x[0]<=0.5+0.05):
        return lambda x,t: 0.0

def vel_u_DBC(x,flag):
    if flag==boundaryTags['bottom'] or flag==boundaryTags['left'] or flag==boundaryTags['right']:
        return lambda x,t: 0.0

def vel_v_DBC(x,flag):
    if flag==boundaryTags['bottom'] or flag==boundaryTags['left'] or flag==boundaryTags['right']:
        return lambda x,t: 0.0
    elif flag == boundaryTags['top'] and (x[0]>=0.5-0.05 and x[0]<=0.5+0.05):
        return lambda x,t: -1.0

# ADVECTIVE FLUX #
def pressure_AFBC(x,flag):
    if not(flag == boundaryTags['top']):
        return lambda x,t: 0.0

def pressure_increment_AFBC(x,flag):
    if not (flag == boundaryTags['top']):
        return lambda x,t: 0.0
    elif x[0]>=0.5-0.05 and x[0]<=0.5+0.05:
        return lambda x,t: -1.0

# DIFFUSIVE FLUX #
def pressure_increment_DFBC(x,flag):
    if not (flag == boundaryTags['top']):
        return lambda x,t: 0.0
    elif x[0]>=0.5-0.05 and x[0]<=0.5+0.05:
        return lambda x,t: 0.0

#############
# LEVEL SET #
#############
def clsvof_DBC(x,flag):
    if flag == boundaryTags['top']:
        if x[0]>=0.5-0.05 and x[0]<=0.5+0.05:
            return lambda x,t: 0.0 #input water
        else:
            return lambda x,t: 1.0 #in open top area let only air in

def clsvof_AFBC(x,flag):
    if flag != boundaryTags['top']:
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
    'pressure_increment_DBC': pressure_increment_DBC,
    'vel_u_DBC': vel_u_DBC,
    'vel_v_DBC': vel_v_DBC,
    'clsvof_DBC': clsvof_DBC,
    # ADVECTIVE FLUX BCs #
    'pressure_AFBC': pressure_AFBC,
    'pressure_increment_AFBC': pressure_increment_AFBC,
    'vel_u_AFBC': lambda x,flag: None,
    'vel_v_AFBC': lambda x,flag: None,
    'clsvof_AFBC': clsvof_AFBC,
    # DIFFUSIVE FLUX BCs #
    'pressure_increment_DFBC': pressure_increment_DFBC,
    'vel_u_DFBC': lambda x,flag: lambda x,t: 0.0,
    'vel_v_DFBC': lambda x,flag: lambda x,t: 0.0,
    'clsvof_DFBC': lambda x,flag: None}
myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=opts.ns_model,
<<<<<<< HEAD
=======
                                             ls_model=1,
>>>>>>> TwoPhaseFlow
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
                                             boundaryConditions=boundaryConditions)
<<<<<<< HEAD
myTpFlowProblem.physical_parameters['densityA'] = 1800.0
myTpFlowProblem.physical_parameters['viscosityA'] = 500.0/1800.0
myTpFlowProblem.physical_parameters['densityB'] = 1.0
myTpFlowProblem.physical_parameters['viscosityB'] = 2.0E-5/1.0
myTpFlowProblem.physical_parameters['surf_tension_coeff'] = 0.
myTpFlowProblem.rans3p_parameters['ns_forceStrongDirichlet']=True
#myTpFlowProblem.clsvof_parameters['lambdaFact']=1.0
=======
myTpFlowProblem.Parameters.physical['densityA'] = 1800.0
myTpFlowProblem.Parameters.physical['viscosityA'] = 500.0/1800.0
myTpFlowProblem.Parameters.physical['densityB'] = 1.0
myTpFlowProblem.Parameters.physical['viscosityB'] = 2.0E-5/1.0
myTpFlowProblem.Parameters.physical['surf_tension_coeff'] = 0.
myTpFlowProblem.Parameters.physical.gravity = np.array([0., -9.8, 0.])
#myTpFlowProblem.clsvof_parameters['lambdaFact']=1.0

myTpFlowProblem.useBoundaryConditionsModule = False
myTpFlowProblem.Parameters.Models.rans3p.epsFact_viscosity = 3.
myTpFlowProblem.Parameters.Models.rans3p.epsFact_density = 3.
myTpFlowProblem.Parameters.Models.rans3p.ns_shockCapturingFactor = 0.5
myTpFlowProblem.Parameters.Models.rans3p.timeDiscretization = 'vbdf'
myTpFlowProblem.Parameters.Models.rans3p.ns_forceStrongDirichlet = True

myTpFlowProblem.outputStepping.systemStepExact = True

myTpFlowProblem.Parameters.mesh.triangleOptions = "VApq30Dena%8.8f" % (old_div((he ** 2), 2.0),)
>>>>>>> TwoPhaseFlow
