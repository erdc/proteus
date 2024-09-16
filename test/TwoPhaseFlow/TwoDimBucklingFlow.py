"""
Multiphase Flow Test
"""
import numpy as np
from proteus import (Domain, Context,
                     MeshTools as mt)
from proteus.Profiling import logEvent
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
import os

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ('ns_model',1,"ns_model = {rans2p,rans3p}"),
    ("final_time",3.0,"Final time for simulation"),
    ("dt_output",0.1,"Time interval to output solution"),
    ("cfl",0.2,"Desired CFL restriction"),
    ("refinement",16,"level of refinement"),
    ("he",1E-2,"he value"),
    ("ARTIFICIAL_VISCOSITY",3,"artificial viscosity")
    ])

assert opts.ns_model==1, "use ns_model=1 (rans3pf) for this"

# ****************** #
# ***** GAUGES ***** #
# ****************** #
# None

# *************************** #
# ***** DOMAIN AND MESH ***** #
# ****************** #******* #
useHex=False
structured=False
nd=2
he=opts.he #Used if unstructured. Goal: he=0.0025=2.5E-3
refinement=opts.refinement #Used if structured. Goal: 64
L=(1.0,1.0)
if useHex:
    raise NotImplementedError
else:
    boundaries = ['bottom', 'right', 'top', 'left', 'front', 'back']
    boundaryTags = dict([(key, i + 1) for (i, key) in enumerate(boundaries)])
    if structured:
        nnx = 4 * refinement #he~0.0039
        nny = nnx
        nnz = None
        triangleFlag=1
        domain = Domain.RectangularDomain(L)
        domain.boundaryTags = boundaryTags
        he = L[0]/(nnx - 1)
    else:
        nnx = None
        nny = None
        vertices = [[0.0, 0.0],  #0
                    [L[0], 0.0],  #1
                    [L[0], L[1]],  #2
                    [0.0, L[1]]]  #3
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
        regions = [[1.2, 0.6]]
        regionFlags = [1]
        domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                      vertexFlags=vertexFlags,
                                                      segments=segments,
                                                      segmentFlags=segmentFlags,
                                                      regions=regions,
                                                      regionFlags=regionFlags)
        #go ahead and add a boundary tags member
        domain.boundaryTags = boundaryTags
        domain.polyfile=os.path.dirname(os.path.abspath(__file__))+"/"+"meshBucklingFlow"
        #domain.writePoly("meshBucklingFlow")
        #domain.MeshOptions.triangleOptions = "VApq30Dena%8.8f" % ((he**2)/2.0,)
        domain.MeshOptions.triangleOptions = "VApq30Dena%8.8f" % ((he ** 2) / 2.0,)
        domain.MeshOptions.he = opts.he
        domain.MeshOptions.triangleFlag=0
        domain.MeshOptions.genMesh = False
        logEvent("""Mesh generated using: tetgen -%s %s""" % (domain.MeshOptions.triangleOptions, domain.polyfile + ".poly"))

# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
class zero(object):
    def uOfXT(self,x,t):
        return 0.

flat_inlet=True
class clsvof_init_cond(object):
    def uOfXT(self,X,t):
        x=X[0]
        y=X[1]
        if flat_inlet:
            if x-L[0]/2. >= -0.05 and x-L[0]/2. <= 0.05:
                if y==L[1]:
                    return 0.
                else:
                    return L[1] - y
            else:
                return min(np.sqrt( (x-(L[0]/2.-0.05))**2 + (y-L[1])**2),
                           np.sqrt( (x-(L[0]/2.+0.05))**2 + (y-L[1])**2) )
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
    if flag == boundaryTags['top'] and not (x[0]>=0.5-0.05 and x[0]<=0.5+0.05):
        return lambda x,t: 0.0

def vel_u_DBC(x,flag):
    if flag==boundaryTags['bottom'] or flag==boundaryTags['left'] or flag==boundaryTags['right']:
        return lambda x,t: 0.0
    elif flag == boundaryTags['top'] and (x[0]>=0.5-0.05 and x[0]<=0.5+0.05):
        return lambda x,t: 0.0

def vel_v_DBC(x,flag):
    if flag==boundaryTags['bottom'] or flag==boundaryTags['left'] or flag==boundaryTags['right']:
        return lambda x,t: 0.0
    elif flag == boundaryTags['top'] and (x[0]>=0.5-0.05 and x[0]<=0.5+0.05):
        return lambda x,t: -1.0
    #return lambda x,t: 100*(x[0]-0.45)*(x[0]-0.55)*1/(0.55*0.45)

# ADVECTIVE FLUX #
def pressure_AFBC(x,flag):
    if not(flag == boundaryTags['top'] and not (x[0]>=0.5-0.05 and x[0]<=0.5+0.05)):
        return lambda x,t: 0.0

def pressure_increment_AFBC(x,flag):
    if not (flag == boundaryTags['top']):
        return lambda x,t: 0.0
    elif x[0]>=0.5-0.05 and x[0]<=0.5+0.05:
        return lambda x,t: -1.0
    #return lambda x,t: 100*(x[0]-0.45)*(x[0]-0.55)*1/(0.55*0.45) # -1.0

# NOTE: recall that D.BCs are set strongly so I want to kill the advective boundary integral
def vel_u_AFBC(x,flag):
    if not (flag == boundaryTags['top'] and not (x[0]>=0.5-0.05 and x[0]<=0.5+0.05)):
        return lambda x,t: 0.0

def vel_v_AFBC(x,flag):
    if not (flag == boundaryTags['top'] and not (x[0]>=0.5-0.05 and x[0]<=0.5+0.05)):
        return lambda x,t: 0.0

# DIFFUSIVE FLUX #
def pressure_increment_DFBC(x,flag):
    if not (flag == boundaryTags['top']):
        return lambda x,t: 0.0
    elif x[0]>=0.5-0.05 and x[0]<=0.5+0.05:
        return lambda x,t: 0.0

vel_u_DFBC = lambda x,flag: lambda x,t: 0.0
vel_v_DFBC = lambda x,flag: lambda x,t: 0.0

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
    if flag == boundaryTags['top']:
        None
    else:
        return lambda x,t: 0.0

clsvof_DFBC = lambda x,flag: None

############################################
# ***** Create myTwoPhaseFlowProblem ***** #
############################################
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
myTpFlowProblem = TpFlow.TwoPhaseFlowProblem()

myTpFlowProblem.domain = domain

myTpFlowProblem.outputStepping.final_time = opts.final_time
myTpFlowProblem.outputStepping.dt_output = opts.dt_output
myTpFlowProblem.SystemPhysics.setDefaults()
myTpFlowProblem.SystemPhysics.useDefaultModels(flowModel=1,interfaceModel=1)

myTpFlowProblem.SystemNumerics.cfl=opts.cfl
myTpFlowProblem.SystemNumerics.useSuperlu=True

myTpFlowProblem.SystemPhysics.modelDict['flow'].p.initialConditions['u']=zero()
myTpFlowProblem.SystemPhysics.modelDict['flow'].p.initialConditions['v']=zero()
myTpFlowProblem.SystemPhysics.modelDict['clsvof'].p.initialConditions['clsvof'] = clsvof_init_cond()
myTpFlowProblem.SystemPhysics.modelDict['pressure'].p.initialConditions['p'] = zero()
myTpFlowProblem.SystemPhysics.modelDict['pressureInc'].p.initialConditions['pInc'] = zero()
myTpFlowProblem.SystemPhysics.modelDict['pressureInit'].p.initialConditions['pInit'] = zero()

myTpFlowProblem.SystemPhysics.boundaryConditions=boundaryConditions

myTpFlowProblem.SystemPhysics.rho_0 = 1800.0
myTpFlowProblem.SystemPhysics.nu_0 = 500.0 /myTpFlowProblem.SystemPhysics.rho_0
myTpFlowProblem.SystemPhysics.rho_1 = 1.0
myTpFlowProblem.SystemPhysics.nu_1 = 2.0E-5/myTpFlowProblem.SystemPhysics.rho_1
myTpFlowProblem.SystemPhysics.surf_tension_coeff = 0.
myTpFlowProblem.SystemPhysics.gravity = np.array([0., -9.8, 0.])

myTpFlowProblem.SystemPhysics.useBoundaryConditionsModule = False
m = myTpFlowProblem.SystemPhysics.modelDict
m['clsvof'].p.coefficients.disc_ICs = False
m['flow'].p.coefficients.useVF=1.0
m['flow'].p.coefficients.forceStrongDirichlet = True
m['flow'].p.coefficients.forceStrongDirichlet = True
m['flow'].p.coefficients.ARTIFICIAL_VISCOSITY = opts.ARTIFICIAL_VISCOSITY
m['flow'].p.coefficients.epsFact_density = 3.
m['flow'].n.ShockCapturingOptions.shockCapturingFactor = 0.5

myTpFlowProblem.outputStepping.systemStepExact = True
