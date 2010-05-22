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
from squareWithHole2dDomain import *
domain = squareWithHole2d(height=1.0,
                          length=1.0,
                          radius=0.25,
                          center=[0.5,0.5],
                          n_points_on_hole=20,
                          cross_section=circular_cross_section)
domain.writePoly("squarehole")
#domain.writePLY("squarehole")
#polyfile = "squarehole"

initialConditions = None

analyticalSolution = None

coefficients = LinearElasticity(E=10.0,nu=0.25,g=[0.0,-9.8],nd=nd)#stokes_2D_tube_p.Stokes2D()
#coefficients = LinearElasticity(E=10.0,nu=0.25,g=[0.0,0.0],nd=nd)#stokes_2D_tube_p.Stokes2D()

def getDBC_u(x,flag):
    if flag == domain.boundaryTags['left']:#x[0] == 0.0:
        return lambda x,t: 0.0
    if flag == domain.boundaryTags['right']:#x[0] == 1.0:
        return lambda x,t: 0.5

def getDBC_v(x,flag):
    if flag == domain.boundaryTags['left']:#x[0] == 0.0:
        return lambda x,t: 0.0
    if flag == domain.boundaryTags['right']:#x[0] == 1.0:
        return lambda x,t: 0.0

dirichletConditions = {0:getDBC_u,
                       1:getDBC_v}

fluxBoundaryConditions = {0:'noFlow',
                          1:'noFlow'}

advectiveFluxBoundaryConditions =  {}

slipBoundaries = [domain.boundaryTags['left'],domain.boundaryTags['left'],domain.boundaryTags['hole']]
def getFluxBC_uu(x,flag):
    if flag in slipBoundaries:
        return lambda x,t: 0.0

def getFluxBC_uv(x,flag):
    if flag in slipBoundaries:
        return lambda x,t: 0.0

def getFluxBC_vu(x,flag):
    if flag in slipBoundaries:
        return lambda x,t: 0.0

def getFluxBC_vv(x,flag):
    if flag in slipBoundaries:
        return lambda x,t: 0.0
    
diffusiveFluxBoundaryConditions = {0:{0:getFluxBC_uu,1:getFluxBC_uv},1:{0:getFluxBC_vu,1:getFluxBC_vv}}

