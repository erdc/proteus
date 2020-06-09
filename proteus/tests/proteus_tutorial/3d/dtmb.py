"""
Multiphase Flow Test
"""
from __future__ import division
from past.utils import old_div
import numpy as np
from proteus import (Domain, Context, Gauges, LinearSolvers,
                     MeshTools as mt)
from proteus.Gauges import PointGauges, LineIntegralGauges, LineGauges
from proteus.Profiling import logEvent
from subprocess import check_call
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus.TwoPhaseFlow.utils.Parameters import ParametersPhysical as PP
from proteus.ctransportCoefficients import smoothedHeaviside
import math

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts= Context.Options([
    ("speed",1.0,"Vessel speed [m/s]"),
    ("final_time",100,"Final time for simulation"),
    ("dt_output",0.5,"Time interval to output solution"),
    ("cfl",0.9,"Desired CFL restriction"),
    ("he",0.3,"Max mesh element diameter"),
    ("skip_gmsh",False,"Assume mesh has already been generated"),
    ])
L = [10, 4, 2]
x0 = [-2, -2, -1]
he = opts.he
boundaryTags={'left':3,'right':4,'bottom':2,'top':1,'front':5,'back':6,'hull':7}
bt = boundaryTags
domain = Domain.PiecewiseLinearComplexDomain()
domain.MeshOptions.setParallelPartitioningType('node')
domain.boundaryTags = boundaryTags
domain.polyfile="dtmb"
domain.geofile="dtmb"
from proteus import Comm
comm = Comm.get()
if not opts.skip_gmsh and comm.isMaster():
    #domain.MeshOptions.triangleOptions="VApq1.25q12feena%e" % ((he**3)/6.0,)
    #gmsh_cmd = "gmsh {0:s} -v 10 -3 -o {1:s} -format msh2".format("boolean.geo", domain.geofile+".msh")
    gmsh_cmd = "gmsh {0:s} -v 10 -3 -o {1:s} -format msh2 -clmax {2:e}".format(domain.geofile+".geo", domain.geofile+".msh",he)
    check_call(gmsh_cmd, shell=True)
    mt.msh2simplex(domain.geofile, nd=3)
    check_call("tetgen -Vfeen {0:s}.ele".format(domain.polyfile), shell=True)
    check_call("mv {0:s}.1.ele {0:s}.ele".format(domain.polyfile), shell=True)
    check_call("mv {0:s}.1.node {0:s}.node".format(domain.polyfile), shell=True)
    check_call("mv {0:s}.1.face {0:s}.face".format(domain.polyfile), shell=True)
    check_call("mv {0:s}.1.neigh {0:s}.neigh".format(domain.polyfile), shell=True)
    check_call("mv {0:s}.1.edge {0:s}.edge".format(domain.polyfile), shell=True)
domain.MeshOptions.genMesh=False

# ****************************** #
# ***** INITIAL CONDITIONS ***** #
# ****************************** #
class zero(object):
    def uOfXT(self,x,t):
        return 0.

# Initial condition
waterLine_x = 1000.0
waterLine_z = 0.25

def signedDistance(x):
    from math import sqrt
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
            return sqrt(phi_x**2 + phi_z**2)

pp = PP()
rho_1 = pp.densityB
rho_0 = pp.densityA
g = [0., 0., -9.81]
class P:
    def __init__(self,waterLevel):
        self.waterLevel=waterLevel
    def uOfXT(self,x,t):
        if signedDistance(x) < 0:
            return -(L[2] - self.waterLevel)*rho_1*g[2] - (self.waterLevel - x[2])*rho_0*g[2]
        else:
            return -(L[2] - x[2])*rho_1*g[2]

class V:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Phi:
    def uOfXT(self,x,t):
        return signedDistance(x)

class VOF:
    def uOfXT(self,x,t):
        return smoothedHeaviside(1.5*he,signedDistance(x))

initialConditions = {'pressure': P(waterLine_z),
                     'vel_u': zero(),
                     'vel_v': zero(),
                     'vel_w': zero(),
                     'vof': VOF(),
                     'rdls': Phi(),
                     'ncls': Phi()}


# ******************************* #
# ***** BOUNDARY CONDITIONS ***** #
# ******************************* #
non_slip_BCs=True
openTop=True
model_scale=24.830
def velRamp(t):
    return (opts.speed/model_scale)*min(1.0,t/40.0)

# DIRICHLET BOUNDARY CONDITIONS #
def vel_u_DBC(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: velRamp(t)*(1.0-initialConditions['vof'].uOfXT(x,t))
    elif flag == boundaryTags['hull']:
        return lambda x,t: 0.0
    elif flag == boundaryTags['right']:
        return lambda x,t: 0.0

def vel_v_DBC(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: 0.0
    elif flag == boundaryTags['hull']:
        return lambda x,t: 0.0
    elif flag == boundaryTags['right']:
        return lambda x,t: 0.0

def vel_w_DBC(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t: 0.0
    elif flag == boundaryTags['hull']:
        return lambda x,t: 0.0
    elif flag == boundaryTags['right']:
        return lambda x,t: 0.0

def pressure_DBC(x,flag):
    if flag == boundaryTags['right']:
        return lambda x,t: initialConditions['pressure'].uOfXT(x,t)

def vof_DBC(x,flag):
    if flag in [boundaryTags['left'], boundaryTags['right']]:
        return lambda x,t: initialConditions['vof'].uOfXT(x,t)

def phi_DBC(x,flag):
    if flag in [boundaryTags['left'], boundaryTags['right']]:
        return lambda x,t: initialConditions['ncls'].uOfXT(x,t)


# ADVECTIVE FLUX BOUNDARY CONDITIONS #
def vel_u_AFBC(x,flag):
    if flag in [0,boundaryTags['front'],boundaryTags['back'],boundaryTags['bottom'],boundaryTags['top']]:
        return lambda x,t: 0.0

def vel_v_AFBC(x,flag):
    if flag in [0,boundaryTags['front'],boundaryTags['back'],boundaryTags['bottom'],boundaryTags['top']]:
        return lambda x,t: 0.0

def vel_w_AFBC(x,flag):
    if flag in [0,boundaryTags['front'],boundaryTags['back'],boundaryTags['bottom'],boundaryTags['top']]:
        return lambda x,t: 0.0

def pressure_AFBC(x,flag):
    if flag  == boundaryTags['left']:
        return lambda x,t: -velRamp(t)*(1.0-initialConditions['vof'].uOfXT(x,t))
    if flag  in [0,boundaryTags['front'],boundaryTags['back'],boundaryTags['bottom'],boundaryTags['top']]:
        return lambda x,t: 0.0

def vof_AFBC(x,flag):
    if flag in [boundaryTags['left'],boundaryTags['right']]:
        return None
    else:
        return lambda x,t: 0.0

# DIFFUSIVE FLUX BOUNDARY CONDITIONS #
def vel_u_DFBC(x,flag):
    if flag in [0,boundaryTags['front'],boundaryTags['back'],boundaryTags['bottom'],boundaryTags['right'],boundaryTags['top']]:
        return lambda x,t: 0.0

def vel_v_DFBC(x,flag):
    if flag in [0,boundaryTags['front'],boundaryTags['back'],boundaryTags['bottom'],boundaryTags['top']]:
        return lambda x,t: 0.0

def vel_w_DFBC(x,flag):
    if flag in [0,boundaryTags['front'],boundaryTags['back'],boundaryTags['bottom'],boundaryTags['top']]:
        return lambda x,t: 0.0

############################################
# ***** Create myTwoPhaseFlowProblem ***** #
############################################
outputStepping = TpFlow.OutputStepping(opts.final_time,dt_output=opts.dt_output)

boundaryConditions = {
    # DIRICHLET BCs #
    'pressure_DBC': pressure_DBC,
    'vel_u_DBC': vel_u_DBC,
    'vel_v_DBC': vel_v_DBC,
    'vel_w_DBC': vel_w_DBC,
    'vof_DBC':  vof_DBC,
    'ncls_DBC': phi_DBC,
    # ADVECTIVE FLUX BCs #
    'pressure_AFBC': pressure_AFBC,
    'vel_u_AFBC': vel_u_AFBC,
    'vel_v_AFBC': vel_v_AFBC,
    'vel_w_AFBC': vel_w_AFBC,
    'vof_AFBC': vof_AFBC,
    # DIFFUSIVE FLUX BCs #
    'vel_u_DFBC': vel_u_DFBC,
    'vel_v_DFBC': vel_v_DFBC,
    'vel_w_DFBC': vel_w_DFBC,
    'vof_DFBC': lambda x, flag: None}

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=0,
                                             ls_model=0,
                                             nd=3,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             he=he,
                                             domain=domain,
                                             initialConditions=initialConditions,
                                             boundaryConditions=boundaryConditions,
                                             useSuperlu=False)
myTpFlowProblem.Parameters.physical.gravity = g
#myTpFlowProblem.Parameters.Models.rans2p.n.linearSmoother = LinearSolvers.NavierStokes_TwoPhasePCD
