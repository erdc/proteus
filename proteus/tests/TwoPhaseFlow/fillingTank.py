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
    ("test_case",1,"1: Guermond's, 2: physical"),
    ('ns_model',1,"ns_model = {rans2p,rans3p}"),
    ("final_time",5.0,"Final time for simulation"),
    ("dt_output",0.2,"Time interval to output solution"),
    ("cfl",0.2,"Desired CFL restriction"),
    ("refinement",3,"level of refinement"),
    ("he",1E-2,"he value"),
    ("ARTIFICIAL_VISCOSITY",3,"artificial viscosity")    
    ])

# ****************** #
# ***** GAUGES ***** #
# ****************** #
# None

# *************************** #
# ***** DOMAIN AND MESH ***** #
# ****************** #******* #
tank_dim = (0.4,0.4) 
# MESH 
refinement = opts.refinement
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
    nnx=nny=nnx=None
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
    he = opts.he
    domain.MeshOptions.he = he
    domain.MeshOptions.triangleOptions = "VApq30Dena%8.8f" % (old_div((he ** 2), 2.0),)

# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
class zero(object):
    def uOfXT(self,x,t):
        return 0.

IC_type=1
class clsvof_init_cond(object):
    def uOfXT(self,x,t):
        if IC_type==0:
            if x[0]==0 and x[1]>=0.3 and x[1]<=0.35:
                return 0.0
            else:
                if x[1]>=0.3 and x[1]<=0.35:
                    return x[0]
                elif x[1]>=0.35:
                    return np.sqrt(x[0]**2 + (x[1]-(0.35))**2)
                else:
                    return np.sqrt(x[0]**2 + (x[1]-(0.3))**2)        
        #
        else:
            if x[0]<0.01 and x[1]>=0.3 and x[1]<=0.35:
                return -1.0
            elif ((x[0]==0.01 and x[1]>=0.3 and x[1]<=0.35) or (x[0]<=0.01 and (x[1]==0.3 or x[1]==0.35))):
                return 0.0
            else:
                return 1.0

# ******************************* #
# ***** BOUNDARY CONDITIONS ***** #
# ******************************* #
#################
# NAVIER STOKES #
#################
inflow_velocity=0.25
eps=0 #opts.he
# DIRICHLET BCs #
def pressure_DBC(x,flag):
    if flag == boundaryTags['top'] and x[0]>=0.3-eps and x[0]<=0.35+eps:
        return lambda x,t: 0.0
    
def pressure_increment_DBC(x,flag):
    if flag == boundaryTags['top'] and x[0]>=0.3-eps and x[0]<=0.35+eps:
        return lambda x,t: 0.0    
    
def vel_u_DBC(x,flag):
    if flag == boundaryTags['bottom'] or flag == boundaryTags['right']: #non-slip at bottom and right
        return lambda x,t: 0.0
    elif flag == boundaryTags['top']:
        if not (x[0]>=0.3-eps and x[0]<=0.35+eps): # non-slip in closed top area
            return lambda x,t: 0.0
    #################
    # LEFT BOUNDARY #
    #################
    elif flag == boundaryTags['left']:
        if x[1]>=0.3-eps and x[1]<=0.35+eps:
            return lambda x,t: inflow_velocity
        #elif x[1]>=0.25: # Slip
        #    None
        else: # Non-slip
            return lambda x,t: 0.0    

def vel_v_DBC(x,flag):
    if flag == boundaryTags['bottom'] or flag == boundaryTags['right']: #non-slip at bottom and right
        return lambda x,t: 0.0
    elif flag == boundaryTags['top']:
        if not (x[0]>=0.3-eps and x[0]<=0.35+eps): # non-slip in closed top area
            return lambda x,t: 0.0
    #################
    # LEFT BOUNDARY #
    #################
    elif flag == boundaryTags['left']:
        if x[1]>=0.3-eps and x[1]<=0.35+eps:
            return lambda x,t: 0.0
        #elif x[1]>=0.25: # Slip
        #    None
        else: # Non-slip
            return lambda x,t: 0.0

# ADVECTIVE FLUX #
def pressure_AFBC(x,flag):
    if not(flag == boundaryTags['top'] and x[0]>=0.3-eps and x[0]<=0.35+eps):
        return lambda x,t: 0.0

def pressure_increment_AFBC(x,flag):
    if flag == boundaryTags['left'] and x[1]>=0.3-eps and x[1]<=0.35+eps:
        return lambda x,t: -inflow_velocity
    elif not (flag == boundaryTags['top'] and x[0]>=0.3-eps and x[0]<=0.35+eps):
        return lambda x,t: 0.0    

# NOTE: recall that D.BCs are set strongly so I want to kill the advective boundary integral
vel_u_AFBC = lambda x,flag: lambda x,t: 0.0
vel_v_AFBC = lambda x,flag: lambda x,t: 0.0
        
# DIFFUSIVE FLUX #
def pressure_increment_DFBC(x,flag):
    if not (flag == boundaryTags['top'] and x[0]>=0.3-eps and x[0]<=0.35+eps):
        return lambda x,t: 0.0

vel_u_DFBC = lambda x,flag: lambda x,t: 0.0
vel_v_DFBC = lambda x,flag: lambda x,t: 0.0

#############
# LEVEL SET #
#############
eps2=opts.he
def clsvof_DBC(x,flag):
    if flag == boundaryTags['left']:
        if x[1]>=0.30-eps and x[1]<=0.35+eps:
            if x[1]>0.30+eps2 and x[1]<0.35-eps2:
                return lambda x,t: 0.0
            else:
                return lambda x,t: 1.0
    elif flag == boundaryTags['top'] and x[0]>=0.3-eps and x[0]<=0.35+eps:
        return lambda x,t: 1.0 # In open area let only air in            

def clsvof_AFBC(x,flag):
    if flag == boundaryTags['top'] and x[0]>=0.3-eps and x[0]<=0.35+eps:
        None
    elif flag == boundaryTags['left'] and x[1]>=0.3-eps and x[1]<=0.35+eps:
        None
    else:
        return lambda x,t: 0.0

clsvof_DFBC = lambda x,flag: None

############################################
# ***** Create myTwoPhaseFlowProblem ***** #
############################################
outputStepping = TpFlow.OutputStepping(opts.final_time,dt_output=opts.dt_output)
outputStepping.systemStepExact = True
initialConditions = {'pressure': zero(),
                     'pressure_increment': zero(),
                     'vel_u': zero(),
                     'vel_v': zero(),
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
    'vel_u_AFBC': vel_u_AFBC,
    'vel_v_AFBC': vel_v_AFBC,
    'clsvof_AFBC': clsvof_AFBC,
    # DIFFUSIVE FLUX BCs #
    'pressure_increment_DFBC': pressure_increment_DFBC,
    'vel_u_DFBC': vel_u_DFBC,
    'vel_v_DFBC': vel_v_DFBC,
    'clsvof_DFBC': clsvof_DFBC}
myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=opts.ns_model,
                                             ls_model=1,
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
                                             useSuperlu=True)
m = myTpFlowProblem.Parameters.Models
m.clsvof.p.CoefficientsOptions['disc_ICs']=False if IC_type==0 else True
m.rans3p.p.CoefficientsOptions['forceStrongDirichlet']=True
m.rans3p.p.CoefficientsOptions['ARTIFICIAL_VISCOSITY']=opts.ARTIFICIAL_VISCOSITY
if opts.test_case==1:
    physical_parameters=myTpFlowProblem.Parameters.physical
    physical_parameters['densityA'] = 1000.0
    physical_parameters['kinematicViscosityA'] = 1.0/physical_parameters['densityA']
    physical_parameters['densityB'] = 1.0
    physical_parameters['kinematicViscosityB'] = 1.8E-2/physical_parameters['densityB']
    physical_parameters['surf_tension_coeff'] = 0.
    physical_parameters['gravity'] = [0.0, -1.0, 0.0]

myTpFlowProblem.Parameters.mesh.triangleOptions = "VApq30Dena%8.8f" % (old_div((he ** 2), 2.0),)
