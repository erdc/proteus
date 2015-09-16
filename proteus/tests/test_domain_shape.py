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

for bc in caisson3D.bc:
    bc.setNoSlip()

for bc in tank.bc:
    bc.setNoSlip()    
tank.bc_top.setOpen()  # BC can be set individually

def newBC(x, t):
    return 0.0

tank.bc_bottom.DBC_u = newBC  # customized functions can still be passed

