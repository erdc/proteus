"""
tilting-flume
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
from proteus.ctransportCoefficients import smoothedHeaviside, smoothedHeaviside_integral

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ("final_time",5.0,"Final time for simulation"),
    ("dt_output",0.01,"Time interval to output solution"),
    ("cfl",0.25,"Desired CFL restriction"),
    ("waterLevel", 0.9,"Height of water column in m"),
    ("tailwater", 0.9,"Height of water column in m"),
    ("g",(0,-9.81,0), "Gravity vector in m/s^2"),
    ("he",0.1, "element diameter"),
    # run time options
    ("dt_fixed", 0.01, "Fixed time step in s"),
    ("dt_init", 0.001 ,"Maximum initial time step in s"),
    ("flowrate",0.5,"unit flowrate for 2d (assuming thickness of 1m)"), 
    ("slope",0.01,"slope of tilting flume"),
    ("flume_length",7.,"length of flume downstream of flow constriction in m"),
    ("flume_height",.7,"height of flume in m")
    ])

he=opts.he
slope=opts.slope
waterLevel=opts.waterLevel
tailwater=opts.tailwater
qu=opts.flowrate

# ----- PHYSICAL PROPERTIES ----- #                                                                 
# Water                                                                                             
rho_0 = 998.2
nu_0 = 1.004e-6
# Air                                                                                               
rho_1 = 1.205
nu_1 = 1.500e-5
# Surface Tension                                                                                   
sigma_01 = 0.0
# Gravity                                                                                           
g = opts.g

# *************************** #
# ***** DOMAIN AND MESH ***** #
# ****************** #******* #

#Function to tilt vertices:
def tilt(vertices, slope):
    verts=vertices
    global theta
    theta = np.arctan(slope)
    axis = 1728 #dimensionless                                                                      
    axism = axis*2.86512/192.81 #in meters                                                          
    for i in range(len(verts)):
        x1 = verts[i][0]
        y1 = verts[i][1]
        distance = axism - x1
        dx = distance - distance*np.cos(theta)
        dx2 = y1*np.sin(theta)
        dy = distance*np.sin(theta)
        dy2 = y1*np.cos(theta)
        x2 = x1 + dx + dx2
        y2 = dy + dy2
        if i == 0:
            tiltedvertices = [[x2,y2]]
        else:
            tiltedvertices.append([x2,y2])
    return tiltedvertices


##Upstream basin                                                                                    
straight = 0.28617+0.4572+2.12175
straightpoints = 40
spacing = 1000*straight/straightpoints #mm
vertices = [[-(0.28617+0.4572),0]]
for i in range(straightpoints):
    newx = -(0.28617+0.4572)+spacing*(i+1)/1000
    vertices.append([newx,0])

##Concave up portion of flow acceleration                                                          
numvertices = 20
LX1 = 1.400
straightend = vertices[len(vertices)-1][0]
for i in range(numvertices):
    newx = straightend+(i+1)*(LX1/(numvertices))
    a = 1.1507*(10**-7)
    newz = a*(((newx-straightend)*1000)**3)/1000
    vertices.append([newx,newz])

##Concave down portion of flow acceleration                                                        
for i in range(numvertices):
    newx = straightend+LX1+(i+1)*(LX1/(numvertices))
    newz = 0.63152-a*(((2.8-(newx-straightend))*1000)**3)/1000
    vertices.append([newx,newz])

##Test section                                                                                     
begindsstraight = vertices[len(vertices)-1][0]
endx = begindsstraight + opts.flume_length
spacing =1000*(endx - begindsstraight)/straightpoints #mm
for i in range(straightpoints):
    newx = begindsstraight+spacing*(i+1)/1000
    vertices.append([newx,0.63152])

##Top Boundary
vertices.append([endx,0.63152+opts.flume_height])
vertices.append([-(0.28617+0.4572),0.63152+opts.flume_height]) 

#Apply Tilt:
vertices = tilt(vertices, slope)

#VertexFlags                                                                                                    
boundaries=['left','outflow','bottom','top']
boundaryTags=dict([[key,i+1] for [i,key] in enumerate(boundaries)])
vertexFlags=[boundaryTags['left']]
for i in range(len(vertices)-4):
    vertexFlags.append(boundaryTags['bottom'])
vertexFlags.append(boundaryTags['outflow'])
vertexFlags.append(boundaryTags['top'])
vertexFlags.append(boundaryTags['top'])

#Segments
segments = [[0,1]]
for i in range(len(vertices)-2):
    segments.append([i+1,i+2])
segments.append((len(vertices)-1,0))

#SegmentFlags
segmentFlags=[boundaryTags['bottom']]
for i in range(len(segments)-4):
    segmentFlags.append(boundaryTags['bottom'])
segmentFlags.append(boundaryTags['outflow'])
segmentFlags.append(boundaryTags['top'])
segmentFlags.append(boundaryTags['left'])

#RegionFlags
regions=[[0.5,0.5]]
regionFlags=[1]

domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                              vertexFlags=vertexFlags,
                                              segments=segments,
                                              segmentFlags=segmentFlags,
                                              regions = regions,
                                              regionFlags = regionFlags,)

domain.MeshOptions.setParallelPartitioningType('node')
domain.boundaryTags = boundaryTags                                                                            
domain.writePoly("mesh")
domain.writePLY("mesh")
domain.writeAsymptote("mesh")
triangleOptions = "VApq30Dena%8.8f" % ((he**2)/2.0,)                                               
logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))                                                                  
domain.MeshOptions.triangleOptions=triangleOptions
bottom=min(vertices, key=lambda x: x[1])[1]
topy=max(vertices, key=lambda x: x[1])[1]

# ****************************** #                                                                  
# ***** User Defined Functions **#                                                                  
# ****************************** #

def signedDistance(X):
    x=X[0]
    y=X[1]
    return y-waterLevel

u = qu*np.sin(theta)
v = qu*np.cos(theta)

epsFact_consrv_heaviside=3.0    
def twpflowVelocity_u(X,t):
    waterspeed = u
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,X[1]-waterLevel)
    return (1.0-H)*waterspeed

def twpflowVelocity_u_flux(X,t):
    waterspeed = u
    H = smoothedHeaviside(epsFact_consrv_heaviside*he,X[1]-waterLevel)
    return -1*(1.0-H)*waterspeed

def intwpflowPressure(X,t):
    p_top = 0.0
    phi_top = topy - waterLevel
    phi = X[1] - waterLevel
    return p_top - g[1]*(rho_0*(phi_top - phi) + \
                         (rho_1 -rho_0) * \
                         (smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_top)
                          -
                          smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))

def outtwpflowPressure(X,t):
    p_top = 0.0
    phi_top = topy - tailwater
    phi = X[1] - tailwater
    return p_top - g[1]*(rho_0*(phi_top - phi) +
                         (rho_1 -rho_0) * \
                         (smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_top)
                          -
                          smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))

def inbcVF(X,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he, X[1] - waterLevel)

def inbcPhi(X,t):
    return X[1] - waterLevel

def outbcVF(X,t):
    return smoothedHeaviside(epsFact_consrv_heaviside*he, X[1] - tailwater)

def outbcPhi(X,t):
    return X[1] - tailwater

#Apply tilt to inflow from bottom                                                                       
c = tilt([[0,0]], slope)
def inflow_(X):
    r = 0.381
    if ((X[0]-c[0][0])**2+(X[1]-c[0][1])**2)**0.5 < r:
        return 1
    else:
        return 0

# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #

class zero(object):
    def uOfXT(self,x,t):
        return 0.0

class VF_IC:
    def uOfXT(self, x, t):
        return smoothedHeaviside(epsFact_consrv_heaviside*he,signedDistance(x))

class PHI_IC:
    def uOfXT(self, x, t):
        return signedDistance(x)
        
# ******************************* #                                                                 
# ***** BOUNDARY CONDITIONS ***** #                                                                 
# ******************************* #                                                                 

non_slip_BCs=True
openTop=True
# DIRICHLET BOUNDARY CONDITIONS #                                                                   
def vel_u_DBC(x,flag):
    if flag != boundaryTags['outflow']:
        if (flag == boundaryTags['bottom'] and inflow_(x)==1):
            return lambda x,t: u
        else:
            return lambda x,t: 0.0
            
def vel_v_DBC(x,flag):
    if flag != boundaryTags['top']:
        if (flag == boundaryTags['bottom'] and inflow_(x)==1):
            return lambda x,t: v
        else:
            return lambda x,t: 0.0

def vel_w_DBC(x,flag):
    return None

def pressure_DBC(x,flag):
        if flag == boundaryTags['top']:
            return lambda x,t: 0.0 
        if flag == boundaryTags['outflow']:
            return outtwpflowPressure

def vof_DBC(x,flag):
    if flag ==boundaryTags['top']:
        return lambda x,t: 1.0
    elif flag == boundaryTags['outflow']:
        return outbcVF
    else:
        return lambda x,t: 0.0

def ncls_DBC(x,flag):
    if flag == boundaryTags['outflow']:
        return outbcPhi
    else:
        return None
    
def rdls_DBC(x, flag):
    pass

# ADVECTIVE FLUX BOUNDARY CONDITIONS #                                                              
def vel_u_AFBC(x,flag):
    return None
def vel_v_AFBC(x,flag):
    return None
def vel_w_AFBC(x,flag):
    return None

def pressure_AFBC(x,flag):
    if flag != boundaryTags['top']:
        if flag != boundaryTags['outflow']:
            if (flag == boundaryTags['bottom'] and inflow_(x)==1):
                return lambda x,t: -qu
            else:
                return lambda x,t: 0.0

def vof_AFBC(x,flag):
    if flag == boundaryTags['top']:
        return None
    elif flag == boundaryTags['outflow']:
        return None
    else:
        return lambda x,t: 0.0

def ncls_AFBC(x,flag):
    return None

def rdls_AFBC(x,flag):
    return None


############################################
# ***** Create myTwoPhaseFlowProblem ***** #
############################################
outputStepping = TpFlow.OutputStepping(opts.final_time,dt_output=opts.dt_output)
initialConditions = {'pressure': zero(),
                     'vel_u': zero(),
                     'vel_v': zero(),
                     'vof': VF_IC(),
                     'ncls': PHI_IC(),
                     'rdls': PHI_IC()}

boundaryConditions = {
    # DIRICHLET BCs #                                                                               
    'pressure_DBC': pressure_DBC,
    'vel_u_DBC': vel_u_DBC,
    'vel_v_DBC': vel_v_DBC,
    'vof_DBC': vof_DBC,
    'ncls_DBC': ncls_DBC,
    'rdls_DBC':rdls_DBC,
    # ADVECTIVE FLUX BCs #                                                                          
    'pressure_AFBC': pressure_AFBC,
    'vel_u_AFBC': vel_u_AFBC,
    'vel_v_AFBC': vel_v_AFBC,
    'vof_AFBC': vof_AFBC,
    'ncls_AFBC': ncls_AFBC,
    'rdls_AFBC': rdls_AFBC,
    # DIFFUSIVE FLUX BCs #                                                                          
    'vel_u_DFBC': lambda x, flag: lambda x,t: 0.0,
    'vel_v_DFBC': lambda x, flag: lambda x,t: 0.0,
    'vof_DFBC': lambda x, flag: None,
    'ncls_DFBC': lambda x, flag: None,
    'rdls_DFBC': lambda x, flag: None,}

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=0,
                                             ls_model=0,
                                             nd=2,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             he=he,
                                             domain=domain,
                                             initialConditions=initialConditions,
                                             boundaryConditions=boundaryConditions,
                                             useSuperlu=False)
physical_parameters = myTpFlowProblem.physical_parameters
myTpFlowProblem.clsvof_parameters['disc_ICs']=False
