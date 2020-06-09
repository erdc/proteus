"""
dambreak 2-D
"""
from __future__ import division
from past.utils import old_div
import numpy as np
from proteus import (Domain, Context)                     
from proteus.Profiling import logEvent
from proteus.mprans.SpatialTools import Tank2D
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus.Gauges import PointGauges, LineIntegralGauges, LineGauges

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ("final_time",3.0,"Final time for simulation"),
    ("dt_output",0.01,"Time interval to output solution"),
    ("cfl",0.25,"Desired CFL restriction"),
    ("he",0.01,"he relative to Length of domain in x"),
    ("refinement",3,"level of refinement")
    ])


# ****************** #
# ***** GAUGES ***** #
# ****************** #
height_gauges1 = LineGauges(gauges=((("phi",),
                                        (((2.724, 0.0, 0.0),
                                          (2.724, 1.8, 0.0)),
                                         ((2.228, 0.0, 0.0),
                                          (2.228, 1.8, 0.0)),
                                         ((1.732, 0.0, 0.0),
                                          (1.732, 1.8, 0.0)),
                                         ((0.582, 0.0, 0.0),
                                          (0.582, 1.8, 0.0)))),),
                                        fileName="height1.csv")

height_gauges2 = LineGauges(gauges=((("phi",),
                                     (((0.0, 0.0, 0.0),
                                       (0.0, 0.0, -0.01)),
                                      ((0.0, 0.0, 0.0),
                                        (3.22, 0.0, 0.0)))),),
                            fileName="height2.csv")

pressure_gauges = PointGauges(gauges=((('p',),
                                      ((3.22, 0.16, 0.0), #P1                                                               
                                       (3.22, 0.584, 0.0), #P3                                                               
                                       (3.22, 0.12, 0.0))),),
                                       fileName="pressure.csv")


# *************************** #
# ***** DOMAIN AND MESH ***** #
# ****************** #******* #
boundaries=['inflow','outflow','bottom','top']
boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])

vertices=[[0.0,0.0],#0                                                
          [3.22,0.],#1                                          
          [3.22,1.8],#2
          [0., 1.8]]#3


vertexFlags=[boundaryTags['inflow'],
             boundaryTags['bottom'],
             boundaryTags['outflow'],
             boundaryTags['top']]

segments=[[0,1],
          [1,2],
          [2,3],
          [3,0]]



segmentFlags=[boundaryTags['bottom'],
              boundaryTags['outflow'],
              boundaryTags['top'],
              boundaryTags['inflow']]


regx=(vertices[1][0]+vertices[0][0])/2
regy=(vertices[-1][1]+vertices[0][1])/2
regions=[[regx, regy]]
regionFlags=[1]

domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                              vertexFlags=vertexFlags,
                                              segments=segments,
                                              segmentFlags=segmentFlags,
                                              regions = regions,
                                              regionFlags = regionFlags,)
#                                              regionConstraints=regionConstraints)
#                                             holes=holes)                                        

he = opts.he                                                                                       

domain.MeshOptions.setParallelPartitioningType('node')
domain.boundaryTags = boundaryTags
#domain.readPoly("mesh")                                                                            
domain.writePoly("mesh")
domain.writePLY("mesh")
domain.writeAsymptote("mesh")
triangleOptions = "VApq30Dena%8.8f" % ((he**2)/2.0,)                                               
logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))
domain.MeshOptions.triangleOptions=triangleOptions



# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
class zero(object):
    def uOfXT(self,x,t):
        return 0.0

waterLine_y = 0.6
waterLine_x = 1.2

class VF_IC:
    def uOfXT(self, x, t):
        if x[0] < waterLine_x and x[1] < waterLine_y:
            return 0.0
        else:
            return 1.0

class PHI_IC:
    def uOfXT(self, x, t):
        phi_x = x[0] - waterLine_x
        phi_y = x[1] - waterLine_y
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


# ******************************* #                                                                  
# ***** BOUNDARY CONDITIONS ***** #                                                                  
# ******************************* #                                                                  
non_slip_BCs=True
openTop=True
# DIRICHLET BOUNDARY CONDITIONS #                                                                    
def vel_u_DBC(x,flag):
    return None
def vel_v_DBC(x,flag):
    return None
def vel_w_DBC(x,flag):
    return None

def pressure_increment_DBC(x,flag):
    if flag == boundaryTags['top'] and openTop:
        return lambda x,t: 0.0

def pressure_DBC(x,flag):
    if flag == boundaryTags['top'] and openTop:
        return lambda x,t: 0.0

def clsvof_DBC(x,flag):
    if openTop and flag == boundaryTags['top']:
        return lambda x,t: 1.0

def vof_DBC(x,flag):
    if openTop and flag == boundaryTags['top']:
        return lambda x,t: 1.0

def ncls_DBC(x,flag):
    return None

def rdls_DBC(x, flag):
    pass


# ADVECTIVE FLUX BOUNDARY CONDITIONS #                                                              
def vel_u_AFBC(x,flag):
    if openTop and flag == boundaryTags['top']:
        return None
    else: #slip                                                              
        return lambda x,t: 0.0

def vel_v_AFBC(x,flag):
    if openTop and flag == boundaryTags['top']:
        return None
    else: #slip                                                              
        return lambda x,t: 0.0

def vel_w_AFBC(x,flag):
    if openTop and flag == boundaryTags['top']:
        return None
    else: #slip                                                             
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

def vof_AFBC(x,flag):
    if openTop and flag == boundaryTags['top']:
        return None
    else:
        return lambda x,t: 0.0

def ncls_AFBC(x,flag):
    return None

def rdls_AFBC(x,flag):
    return None

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
                     'vof': VF_IC(),
                     'ncls': PHI_IC(),
                     'rdls': PHI_IC()}

boundaryConditions = {
    # DIRICHLET BCs #                                                                                
    'pressure_DBC': pressure_DBC,
#    'pressure_increment_DBC':  pressure_increment_DBC,
    'vel_u_DBC': vel_u_DBC,
    'vel_v_DBC': vel_v_DBC,
#    'vel_w_DBC': vel_w_DBC,
    'vof_DBC': vof_DBC,
    'ncls_DBC': ncls_DBC,
    'rdls_DBC':rdls_DBC,
#    'clsvof_DBC': clsvof_DBC,
    # ADVECTIVE FLUX BCs #                                                                           
    'pressure_AFBC': pressure_AFBC,
#    'pressure_increment_AFBC': pressure_increment_AFBC,
    'vel_u_AFBC': vel_u_AFBC,
    'vel_v_AFBC': vel_v_AFBC,
#    'vel_w_AFBC': vel_w_AFBC,
    'vof_AFBC': vof_AFBC,
    'ncls_AFBC': ncls_AFBC,
    'rdls_AFBC': rdls_AFBC,
#    'clsvof_AFBC': clsvof_AFBC,
    # DIFFUSIVE FLUX BCs #                                                                           
#    'pressure_increment_DFBC': pressure_increment_DFBC,
    'vel_u_DFBC': lambda x, flag: lambda x,t: 0.0,
    'vel_v_DFBC': lambda x, flag: lambda x,t: 0.0,
#    'vel_w_DFBC': lambda x, flag: lambda x,t: 0.0,
    'vof_DFBC': lambda x, flag: None,
    'ncls_DFBC': lambda x, flag: None,
    'rdls_DFBC': lambda x, flag: None,
#    'clsvof_DFBC': lambda x, flag: None
}

#auxVariables={'clsvof': [height_gauges1, height_gauges2],
#              'pressure': [pressure_gauges]}

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=0,
                                             ls_model=0,
                                             nd=2,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             he=he,
#                                             nnx=nnx,
#                                             nny=nny,
#                                             nnz=None,
                                             domain=domain,
                                             initialConditions=initialConditions,
                                             boundaryConditions=boundaryConditions,
#                                             auxVariables=auxVariables,
                                             useSuperlu=False)
physical_parameters = myTpFlowProblem.physical_parameters
#myTpFlowProblem.clsvof_parameters['disc_ICs']=False
