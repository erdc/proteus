#transport of gaussian in 2d periodic domain
nd = 2#number of space dimensions
import numpy
#constant velocity
velocity = numpy.array([1.0,1.0])
#where gaussian starts
center = numpy.array([0.25,0.25])#numpy.array([0.25,0.25])#numpy.array([0.25,0.5])
#size
sigma  = 1./16.
#quadrature
space_quad_order = 4#6 for dgp3
#form of level set equation
useHJ=False#True#False
#number nodes in each direction
nn=11#81#161#41
#number of meshes in multilevel mesh
nLevels = 1
#end time of simulation
T = 2.0
#number of output time steps, ignored for adaptive/CFL based runs
nDTout = 100
#max CFL
runCFL = 0.185#0.3,0.185,0.125 for dgp1,dgp2,dgpk(3)

name = 'la_periodicGauss_2d_dgp2_HJnn21'
