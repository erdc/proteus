"""
Testing module for Domain.py, SHape.py, BC.py
Work in progress
"""


from proteus import Domain
from proteus import Shape as sp
from proteus.default_n import *


# ----- DOMAIN ----- #

domain = Domain.PiecewiseLinearComplexDomain()


# ----- SHAPES ----- #

tank_dim = [4., 4., 4.]
tank_coords = [0., 0., 0.]
tank = sp.Cuboid(domain, tank_dim)
tank.setPosition(tank_coords)

caisson_dim = [.3, .9, .1]
caisson_coords = [2., 2., 2.]
caisson3D = sp.Cuboid(domain, caisson_dim)
caisson3D.setPosition(caisson_coords)
caisson3D.setConstraints([0., 0., 1.], [1., 1., 0.])


# ----- BOUNDARY CONDITIONS ----- #

for bc in caisson3D.BC_list:
    bc.setNoSlip()

for bc in tank.BC_list:
    bc.setNoSlip()
tank.BC_dict['top'].setFreeSlip()  # BC can be set individually by using a dictionary
tank.BC.top.setFreeSlip()  # or by using a class variable (same effect)


boundaryTags = { 'bottom': 1, 'front':2, 'right':3, 'back': 4, 'left':5, 'top':6, 'obstacle':7}
vertices=[[0.,0.,0.],#0
          [1.,0.,0.],#1
          [1.,1.,0.],#2
          [0.,1.,0.],#3
          [0.,0.,1.],#4
          [1.,0.,1.],#5
          [1.,1.,1.],#6
          [0.,1.,1.]]#7
vertexFlags=[boundaryTags['left'],
             boundaryTags['right'],
             boundaryTags['right'],
             boundaryTags['left'],
             boundaryTags['left'],
             boundaryTags['right'],
             boundaryTags['right'],
             boundaryTags['left']]
facets=[[[0,1,2,3]],
        [[0,1,5,4]],
        [[1,2,6,5]],
        [[2,3,7,6]],
        [[3,0,4,7]],
        [[4,5,6,7]]]
facetFlags=[boundaryTags['bottom'],
            boundaryTags['front'],
            boundaryTags['right'],
            boundaryTags['back'],
            boundaryTags['left'],
            boundaryTags['top']]
custom = sp.CustomShape(domain=domain,
                        vertices=vertices,
                        vertexFlags=vertexFlags,
                        facets=facets,
                        facetFlags=facetFlags,
                        boundaryTags=boundaryTags)


def newBC(x, t):
    return 0.0

custom.BC.bottom.u_dirichlet = newBC  # customized functions can still be passed

