"""
Rising bubble test
"""
from __future__ import division
from past.utils import old_div
import numpy as np
from proteus import (Domain, Context)
from proteus.Profiling import logEvent
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ('ns_model',1,"ns_model = {rans2p,rans3p}"),
    ("final_time",5.0,"Final time for simulation"),
    ("dt_output",0.1,"Time interval to output solution"),
    ("cfl",0.2,"Desired CFL restriction"),
    ("he", 0.5, "Maximum mesh element diameter"),
    ("ARTIFICIAL_VISCOSITY",3,"artificial viscosity")
    ])

assert opts.ns_model==1, "Surface tension is only implemented with rans3p. use ns_model=1"

# ****************** #
# ***** GAUGES ***** #
# ****************** #
# None

# *************************** #
# ***** DOMAIN AND MESH ***** #
# ****************** #******* #
L      = [3.22,1.0,1.0]
box_L  = [0.161,0.403,0.35] #0.161
box1_xy = [0.6635,0.2985]
box2_xy = [2.3955,0.2985]
he = opts.he
boundaries=['left','right','bottom','top','front','back','box1_left','box1_right','box1_top','box1_front','box1_back','box2_left','box2_right','box2_top','box2_front','box2_back']
boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
bt = boundaryTags
holes = [[0.5*box_L[0]+box1_xy[0],0.5*box_L[1]+box1_xy[1],0.5*box_L[2]],
         [0.5*box_L[0]+box2_xy[0],0.5*box_L[1]+box2_xy[1],0.5*box_L[2]]]
vertices=[[0.0,0.0,0.0],#0
          [L[0],0.0,0.0],#1
          [L[0],L[1],0.0],#2
          [0.0,L[1],0.0],#3
          [0.0,0.0,L[2]],#4
          [L[0],0.0,L[2]],#5
          [L[0],L[1],L[2]],#6
          [0.0,L[1],L[2]],#7
          [box1_xy[0],box1_xy[1],0.0],#8
          [box1_xy[0]+box_L[0],box1_xy[1],0.0],#9
          [box1_xy[0]+box_L[0],box1_xy[1]+box_L[1],0.0],#10
          [box1_xy[0],box1_xy[1]+box_L[1],0.0],#11
          [box1_xy[0],box1_xy[1],box_L[2]],#12
          [box1_xy[0]+box_L[0],box1_xy[1],box_L[2]],#13
          [box1_xy[0]+box_L[0],box1_xy[1]+box_L[1],box_L[2]],#14
          [box1_xy[0],box1_xy[1]+box_L[1],box_L[2]],#15
          [box2_xy[0],box2_xy[1],0.0],#16
          [box2_xy[0]+box_L[0],box2_xy[1],0.0],#17
          [box2_xy[0]+box_L[0],box2_xy[1]+box_L[1],0.0],#18
          [box2_xy[0],box2_xy[1]+box_L[1],0.0],#19
          [box2_xy[0],box2_xy[1],box_L[2]],#20
          [box2_xy[0]+box_L[0],box2_xy[1],box_L[2]],#21
          [box2_xy[0]+box_L[0],box2_xy[1]+box_L[1],box_L[2]],#22
          [box2_xy[0],box2_xy[1]+box_L[1],box_L[2]]]#23
vertexFlags=[boundaryTags['left'],
             boundaryTags['right'],
             boundaryTags['right'],
             boundaryTags['left'],
             boundaryTags['left'],
             boundaryTags['right'],
             boundaryTags['right'],
             boundaryTags['left'],
             boundaryTags['box1_left'],
             boundaryTags['box1_left'],
             boundaryTags['box1_left'],
             boundaryTags['box1_left'],
             boundaryTags['box1_left'],
             boundaryTags['box1_left'],
             boundaryTags['box1_left'],
             boundaryTags['box1_left'],
             boundaryTags['box2_left'],
             boundaryTags['box2_left'],
             boundaryTags['box2_left'],
             boundaryTags['box2_left'],
             boundaryTags['box2_left'],
             boundaryTags['box2_left'],
             boundaryTags['box2_left'],
             boundaryTags['box2_left']]
facets=[[[0,1,2,3],[8,9,10,11],[16,17,18,19]], # bottom face
        [[0,1,5,4]], # front face
        [[1,2,6,5]], # Right
        [[2,3,7,6]], # back
        [[3,0,4,7]], # left
        [[4,5,6,7]], # Top
        [[8,9,13,12]], # front box 1
        [[9,10,14,13]], # right box 1
        [[10,11,15,14]], # back box 1
        [[11,8,12,15]],  # left box 1
        [[12,13,14,15]], # top box 1
        [[16,17,21,20]], #front box 2
        [[17,18,22,21]], #right box 2
        [[18,19,23,22]], #back box 2
        [[19,16,20,23]], #left box 2
        [[20,21,22,23]]] #top box 2
facetFlags=[boundaryTags['bottom'],
            boundaryTags['front'],
            boundaryTags['right'],
            boundaryTags['back'],
            boundaryTags['left'],
            boundaryTags['top'],
            boundaryTags['box1_front'],
            boundaryTags['box1_right'],
            boundaryTags['box1_back'],
            boundaryTags['box1_left'],
            boundaryTags['box1_top'],
            boundaryTags['box2_front'],
            boundaryTags['box2_right'],
            boundaryTags['box2_back'],
            boundaryTags['box2_left'],
            boundaryTags['box2_top']]
regions=[[0.5*L[0],0.5*L[1],0.5*L[2]]]
regionFlags=[0]

domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                             vertexFlags=vertexFlags,
                                             facets=facets,
                                             facetFlags=facetFlags,
                                             regions = regions,
                                             regionFlags = regionFlags,
                                             holes=holes)
#go ahead and add a boundary tags member
domain.MeshOptions.setParallelPartitioningType('node')
domain.boundaryTags = boundaryTags
domain.writePoly("mesh")
domain.writePLY("mesh")
domain.writeAsymptote("mesh")
domain.MeshOptions.triangleOptions = triangleOptions="VApq1.25q12feena%e" % ((he**3)/6.0,)

# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
class zero(object):
    def uOfXT(self,x,t):
        return 0.

class clsvof_init_cond(object):
    def uOfXT(self,x,t):
        waterLine_z = 0.3
        if x[2]-waterLine_z > 0:
            return 1.0
        elif x[2]-waterLine_z < 0:
            return -1.0
        else:
            return 0.0

#################
# NAVIER STOKES #
#################
# DIRICHLET BCs #
def vel_u_DBC(x,flag):
    if flag == boundaryTags['top']:
        None
    elif flag in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]:
        return lambda x,t: 0.0

def vel_v_DBC(x,flag):
    if flag == boundaryTags['top']:
        None
    elif flag in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]:
        return lambda x,t: 0.0

def vel_w_DBC(x,flag):
    if flag == boundaryTags['top']:
        None
    elif flag in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]:
        return lambda x,t: 0.0

# ADVECTIVE FLUX #
# DIFFUSIVE FLUX #
def vel_u_DFBC(x,flag):
    #if flag == boundaryTags['top']:
        return lambda x,t: 0.0

def vel_v_DFBC(x,flag):
    #if flag == boundaryTags['top']:
        return lambda x,t: 0.0

def vel_w_DFBC(x,flag):
    #if flag == boundaryTags['top']:
        return lambda x,t: 0.0

######################
# PRESSURE INCREMENT #
######################
def pressure_increment_DBC(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0

def pressure_increment_AFBC(x,flag):
    if not flag == boundaryTags['top']:
        return lambda x,t: 0.0

def pressure_increment_DFBC(x,flag):
    if not flag == boundaryTags['top']:
        return lambda x,t: 0.0

############
# PRESSURE #
############
def pressure_DBC(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0

def pressure_AFBC(x,flag):
    if not flag == boundaryTags['top']:
        return lambda x,t: 0.0

##########
# CLSVOF #
##########
def clsvof_DBC(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 1.0 # let only air in
#
def clsvof_AFBC(x,flag):
    if flag == boundaryTags['top']:
        return None
    else:
        return lambda x,t: 0.0

###############
# FORCE FIELD #
###############
def forcex(X,t):
    x = X[0]
    yT = X[2]/L[2]

    # PARAMETERS FOR FORCE #
    wF = 0.5 #width of the force
    bottomForce = 10.0
    topForce = 1.0
    forceStrength = 15.0
    timeModulator = 0.5*(1+np.tanh(t))
    #timeModulator = 1

    # COMPUTE FORCE COEFFICIENTS #
    rightForce =  forceStrength * ((1.0-yT)*bottomForce + yT*topForce)
    leftForce = -forceStrength * ((1.0-yT)*bottomForce + yT*topForce)

    # COMPUTE FORCES #
    rF = rightForce*(x - L[0]/2. >= 0.)*(x - L[0]/2. <= wF/2.)
    lF = leftForce*(x -L[0]/2. <= 0.)*(-wF/2. <= x-L[0]/2)

    return timeModulator*(rF + lF)

def forcey(X,t):
    return 0.

def forcez(X,t):
    return 0.

forceTerms = {0:forcex,
              1:forcey,
              2:forcez}

############################################
# ***** Create myTwoPhaseFlowProblem ***** #
############################################
outputStepping = TpFlow.OutputStepping(opts.final_time,dt_output=opts.dt_output)
outputStepping.systemStepExact = True
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
    'vel_u_AFBC': lambda x, flag: None,
    'vel_v_AFBC': lambda x, flag: None,
    'vel_w_AFBC': lambda x, flag: None,
    'clsvof_AFBC': clsvof_AFBC,
    # DIFFUSIVE FLUX BCs #
    'pressure_increment_DFBC': pressure_increment_DFBC,
    'vel_u_DFBC': vel_u_DFBC,
    'vel_v_DFBC': vel_v_DFBC,
    'vel_w_DFBC': vel_w_DFBC,
    'clsvof_DFBC': lambda x, flag: None}
myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=opts.ns_model,
                                             ls_model=1,
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
                                             boundaryConditions=boundaryConditions,
                                             useSuperlu=True)
myTpFlowProblem.useBoundaryConditionsModule = False
myTpFlowProblem.Parameters.physical['gravity'] = [0.0,0.0,-9.8]
m = myTpFlowProblem.Parameters.Models
m.clsvof.p.CoefficientsOptions['disc_ICs']=True
m.rans3p.p.CoefficientsOptions['ARTIFICIAL_VISCOSITY']=opts.ARTIFICIAL_VISCOSITY
m.rans3p.p.CoefficientsOptions.forceTerms = forceTerms
m.rans3p.n.ShockCapturingOptions.shockCapturingFactor = 0.5

myTpFlowProblem.Parameters.mesh.setParallelPartitioningType('node')
myTpFlowProblem.Parameters.mesh.triangleOptions = triangleOptions="VApq1.25q12feena%e" % ((he**3)/6.0,)
