from pyadh import *
from pyadh.default_p import *
"""
A linear elastic cantilever beam problem
"""
##\page Tests Test Problems 
# \ref linear_elasticity_cantilever_beam_3d_p.py "A linear elastic cantilever beam problem"
#

##\ingroup test
#\file linear_elasticity_cantilever_beam_3d_p.py
#\brief A linear elastic cantilever beam problem
#
nd = 3

L=(1.0,1.0,1.0)

initialConditions = None

analyticalSolution = None

coefficients = LinearElasticity(E=10.0,nu=0.25,g=[0.0,0.0,-9.8],nd=nd)

def getDBC_u(x):
    if x[0] == 0.0:
        return lambda x,t: 0.0
#    if x[0] == 1.0:
#        return lambda x,t: 0.0

def getDBC_v(x):
    if x[0] == 0.0:
        return lambda x,t: 0.0
#    if x[0] == 1.0:
#        return lambda x,t: 0.0

def getDBC_w(x):
    if x[0] == 0.0:
        return lambda x,t: 0.0
#    if x[0] == 1.0:
#        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_u,
                       1:getDBC_v,
                       2:getDBC_w}

fluxBoundaryConditions = {0:'noFlow',
                          1:'noFlow',
                          2:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {}

