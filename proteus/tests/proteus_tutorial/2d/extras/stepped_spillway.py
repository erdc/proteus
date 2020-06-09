"""
Stepped spillway
"""
from __future__ import division
from past.utils import old_div
import numpy as np
from proteus import (Domain, Context)                     
from proteus.Profiling import logEvent
from proteus.mprans.SpatialTools import Tank2D
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus.ctransportCoefficients import smoothedHeaviside, smoothedHeaviside_integral
import stepped_spillway_geom as ssg
#from proteus.Gauges import PointGauges, LineIntegralGauges, LineGauges

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts=Context.Options([
    # water column                                                                          
    ("water_level", 457.0,"Height of water column in m"),
    ("tailwater", 435.0,"Height of water column in m"),
    ("water_width",304, "width (along x) of  water column in m"),
    # tank                                                                                  
    ("tank_dim", (3.22, 1.8), "Dimensions of the tank  in m"),
    #gravity                                                                                
    ("g",(0,-9.81,0), "Gravity vector in m/s^2"),
    # gauges                                                                                
    ("gauge_output", True, "Produce gauge data"),
    ("gauge_location_p", (3.22, 0.12, 0), "Pressure gauge location in m"),
    # mesh refinement and timestep                                                          
    ("refinement", 1 ,"Refinement level, he = L/(4*refinement - 1), where L is the horizont\
al dimension"),
    ("he",0.1, "element diameter"),
    ("cfl", 0.33 ,"Target cfl"),
    # run time options                                                                      
    ("T", 10. ,"Simulation time in s"),
    ("dt_fixed", 0.01, "Fixed time step in s"),
    ("dt_init", 0.001 ,"Maximum initial time step in s"),
    ("useHex", False, "Use a hexahedral structured mesh"),
    ("structured", False, "Use a structured triangular mesh"),
    ("waveheight", 1.0, "wave height"),
    ("wavelength", 10.0, "wave length"),
    ("wallheight",5.0, "levee height"),
    ("HtoV",3.0,"levee horizontal to vertical slope H:V"),
    ("leeward_dry",True, "no water on leeward side"),
    ("gen_mesh", True ,"Generate new mesh"),
    ("test",False, "debugging"),
    ("upstream_length", 20, "distance from inflow boundary to weir"),
    ("downstream_length", 20, "distance from plunge pool to outflow boundary"),
    ("top", 471, "height of the top boundary"),
    ("flowrate",1,"unit flowrate for 2d (assuming thickness of 1m)"),
    ("dt_output",0.01,"output timestep"),
    ])

he=opts.he

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

# ----- CONTEXT ------ #                                                                    

# water                                                                                     
tailwater=opts.tailwater

# *************************** #
# ***** DOMAIN AND MESH ***** #
# ****************** #******* #
boundaries=['bottom','outflow','top','inflow']
boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])

vertices=np.genfromtxt("stepped_spillway_domain.csv", delimiter=",").tolist()
vertices=ssg.ups_len(vertices,opts.upstream_length)
vertices=ssg.dwns_len(vertices,opts.downstream_length)
vertices=ssg.top(vertices,opts.top)
regx=(vertices[1][0]+vertices[0][0])/2
regy=(vertices[-1][1]+vertices[0][1])/2

vertexFlags=ssg.v_flags(vertices,boundaryTags)

segments=ssg.segs(vertices)
segmentFlags=ssg.s_flags(segments,boundaryTags)
aa=1#70                                                                                           
bb=315#int(min(range(len(vertices)), key=lambda i: abs(np.array(vertices)[i,1]-tailwater)))       
refverts=[[vertices[aa][0],vertices[aa][1]+4],[vertices[bb][0],vertices[bb][1]+2]]

refx=(vertices[int((aa+bb)*0.5)][0]*1.01)
refy=(vertices[int((aa+bb)*0.5)][1])
refFlags=[0,0,0,0]
vertices=vertices+refverts
vertexFlags=vertexFlags+refFlags

vn=len(vertices)
refSegs=[[aa,vn-2],[vn-2,vn-1],[vn-1,bb]]
segments=segments+refSegs
segmentFlags=segmentFlags+refFlags
regions=[[regx, regy],[refx,refy]]
regionFlags=[1,2]
#he=(he**2)/2.0                                                                                   
regionConstraints=[he,0.1*he]


domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                              vertexFlags=vertexFlags,
                                              segments=segments,
                                              segmentFlags=segmentFlags,
                                              regions = regions,
                                              regionFlags = regionFlags,
                                              regionConstraints=regionConstraints)
   

domain.MeshOptions.setParallelPartitioningType('node')
domain.boundaryTags = boundaryTags
#domain.readPoly("mesh")                                                                           
domain.writePoly("mesh")
domain.writePLY("mesh")
domain.writeAsymptote("mesh")
triangleOptions = "VApq30Dena%8.8f"# % ((he**2)/2.0,)                                              
logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))
#proteus.MeshTools.                                                                                
domain.MeshOptions.triangleOptions=triangleOptions

waterLevel= opts.water_level
bottom=min(vertices, key=lambda x: x[1])[1]
topy=max(vertices, key=lambda x: x[1])[1]

# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #

weir=vertices[40][0]#opts.water_width                                                             
xstep=vertices[int(min(range(len(vertices)), key=lambda i: abs(np.array(vertices)[i,1]-tailwater)))][0]

def signedDistance(X):
    x=X[0]
    y=X[1]
    if x < weir:
        return y-waterLevel
    elif x > xstep:
        return y-tailwater
    else:
        d1= ((y-waterLevel)**2 + (x-weir)**2)**0.5
        d2= ((y-tailwater)**2 + (x-xstep)**2)**0.5
        return min(d1,d2)

class zero(object):
    def uOfXT(self,x,t):
        return 0.0

waterLine_y = 0.6
waterLine_x = 1.2


class PHI_IC:
    def uOfXT(self, x, t):
        return signedDistance(x)

epsFact_consrv_heaviside = 3.0
class VF_IC:
    def uOfXT(self, x, t):
        return smoothedHeaviside(epsFact_consrv_heaviside*he,signedDistance(x))

# ******************************* #
# ***** BC input functions  ***** #
# ******************************* #

def twpflowVelocity_u(X,t):
   waterspeed = u
   H = smoothedHeaviside(epsFact_consrv_heaviside*he,X[1]-waterLevel)
   return (1.0-H)*waterspeed

u=opts.flowrate/(opts.water_level-vertices[0][1])

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
    return p_top - g[1]*(rho_0*(phi_top - phi) + \
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

# ******************************* #                                                                  
# ***** BOUNDARY CONDITIONS ***** #                                                                  
# ******************************* #                                                                  

non_slip_BCs=True
openTop=True
# DIRICHLET BOUNDARY CONDITIONS #                                                                    
def vel_u_DBC(x,flag):
    if flag == boundaryTags['inflow']:
        return twpflowVelocity_u
    if flag != boundaryTags['outflow']:
        return lambda x,t: 0.0

def vel_v_DBC(x,flag):
    if flag != boundaryTags['top']:
        return lambda x,t: 0.0

def vel_w_DBC(x,flag):
    return None

def pressure_increment_DBC(x,flag):
    return None

def pressure_DBC(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0
    if flag == boundaryTags['outflow']:
        return outtwpflowPressure

def clsvof_DBC(x,flag):
    if openTop and flag == boundaryTags['top']:
        return lambda x,t: 1.0

def vof_DBC(x,flag):
    if flag ==boundaryTags['top']:
        return lambda x,t: 1.0
    elif flag == boundaryTags['inflow']:
        return inbcVF
    elif flag == boundaryTags['outflow']:
        return outbcVF
    else:
        return lambda x,t: 1.0

def ncls_DBC(x,flag):
    if flag == boundaryTags['inflow']:
        return inbcPhi
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
            if flag == boundaryTags['inflow']:
                return twpflowVelocity_u_flux
            else:
                return lambda x,t: 0.0

def clsvof_AFBC(x,flag):
    if openTop and flag == boundaryTags['top']:
        return None
    else:
        return lambda x,t: 0.0

def vof_AFBC(x,flag):
    if flag == boundaryTags['top']:
        return None
    elif flag == boundaryTags['inflow']:
        return None
    elif flag == boundaryTags['outflow']:
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
outputStepping = TpFlow.OutputStepping(opts.T,dt_output=opts.dt_output)
initialConditions = {'pressure': zero(),
#                     'pressure_increment': zero(),
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
#    'clsvof_AFBC': clsvof_AFBC,
    # DIFFUSIVE FLUX BCs #                                                                           
    'vel_u_DFBC': lambda x, flag: lambda x,t: 0.0,
    'vel_v_DFBC':  lambda x, flag: lambda x,t: 0.,
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
