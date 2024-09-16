#! /usr/bin/env python

"""
common definitions for simple transient groundwater problem

    
    
    S_s\pd{h}{t} + \deld \vec q = s(x,t)

          \vec q =-\frac{K(x,t)}{\hat{\mu}(c(x,t))} \grad h

Here K is the Saturated Conductivity, K = k\varrho_0 |g|/\mu_0 and  S_s is the specific storage

Note that we assume S_h are independent of h and t and approximate
the accumulation term in the 'conservative form' 

    \pd{m}{t} = \pd{S_h h}{t}

since this is the form that proteus uses by default
     
"""
from proteus import *
from proteus import SubsurfaceTransportCoefficients as STC

import numpy as np


name = "TransientSinglePhaseFlow"

### Spatial Domain ###
L = (1000,1000,1) #m
nd= 2


#temporal domain
#for converting from seconds to days
days2sec = 60.*60.*24.

ndays = 5000.0#3650.0
T = ndays #[d]
nDTout = 10

#just a simple rectangular domain for now
domain = Domain.RectangularDomain(L[:nd],
                                          units="m")
refinement_level = 16 #define characteristic length
he = L[0]/float(refinement_level)

### Material Properties ###
#homogeneous for now
nmaterial = 1
K_0 = 4.0# m/d
S_s0= 1.0e-2
f_0 = 0.0e-2#1/d 

Ident = np.zeros((nd,nd),'d')
Ident[0,0]=1.0; Ident[1,1] = 1.0

#dictionaries that return functions giving material properties as a function of x,t
#the key to the dictionary is the material id
material_id_base = 0
conductivities = {material_id_base:lambda x,t: K_0*Ident}
sources = {material_id_base:lambda x,t: f_0}
S_ss = {material_id_base:lambda x,t: S_s0}


### Boundary Conditions ###
## Flow Equation
#Pressure boundary condition on bottom left and top right
#no flow everywhere else
inflow_length = L[1]*0.1
outflow_length= L[1]*0.8

#piezometric head
head_inflow = 100.0 #[m]
head_outflow= 0.0

#
bndry_eps=1.0e-6
def head_bc(x,flag):
    if x[0] <= bndry_eps and x[1] <= inflow_length:
        return lambda x,t: head_inflow
    if x[0] >= L[0]-bndry_eps and x[1] >= outflow_length:
        return lambda x,t: head_outflow

def noflux(x,flag):
    if x[1] <= bndry_eps or x[1] >= L[1]-bndry_eps:
        return lambda x,t: 0.0
    if x[0] <= bndry_eps and x[1] > inflow_length:
        return lambda x,t: 0.0
    if x[0] >= L[0]-bndry_eps and x[1] < outflow_length:
        return lambda x,t: 0.0


### Initial conditions ###
class ConstantIC(object):
    def __init__(self,val=0.0):
        self.val = val
    def uOfXT(self,x,t):
        return self.val

##Flow Equation
initialConditions_flow = {0:ConstantIC(head_outflow)}

# numerics
parallel = False

domain.MeshOptions.nnx = domain.MeshOptions.nny = int(L[0]/he)
nLevels = 1
if parallel:
    nLevels = 1

# parallel partitioning
if parallel:
    nLayersOfOverlapForParallel = 1     
    parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.node
    numericalFluxType = NumericalFlux.Advection_DiagonalUpwind_Diffusion_SIPG_exterior

# 
gw_quad_order = 3
