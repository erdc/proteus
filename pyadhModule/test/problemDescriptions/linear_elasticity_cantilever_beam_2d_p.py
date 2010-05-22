from pyadh import *
from pyadh.default_p import *

"""
A linear elastic cantilever beam problem
"""
##\page Tests Test Problems 
# \ref linear_elasticity_cantilever_beam_2d_p.py "A linear elastic cantilever beam problem"
#

##\ingroup test
#\file linear_elasticity_cantilever_beam_2d_p.py
#\brief A linear elastic cantilever beam problem
#
nd = 2

#L=(1.0,1.0,1.0)
from box2dDomain import *
domain = box2D(Lx=10.0,Ly=1.0)
domain.writePoly("cantileverBeam")

initialConditions = None

analyticalSolution = None

coefficients = LinearElasticity(E=317986,nu=0.31,g=[0.0,-9.8],nd=nd)#stokes_2D_tube_p.Stokes2D()

bfs = domain.boundaryFlags
def getDBC_u(x,flag):
    if flag == bfs['left']:
        return lambda x,t: 0.0

def getDBC_v(x,flag):
    if flag == bfs['left']:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_u,
                       1:getDBC_v}

fluxBoundaryConditions = {0:'noFlow',
                          1:'noFlow'}

advectiveFluxBoundaryConditions =  {}

def noTraction(x,flag):
    if flag in [bfs['right'],bfs['front'],bfs['back']]:
        return lambda x,t: 0.0

diffusiveFluxBoundaryConditions = {0:{0:noTraction,1:noTraction},
                                   1:{0:noTraction,1:noTraction}}

