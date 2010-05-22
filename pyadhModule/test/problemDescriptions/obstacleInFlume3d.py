"""
A helper module for doing air/water flow around a moving rigid cylinder in 2D
"""
from math import *
from pyadh import *

#
#material properties
#
#water
rho_0=1000.0#998.2
nu_0=0.001#1.004e-6
#air
rho_1=1.0#1.205
nu_1= 0.00002#1.500e-5
#rho_1=rho_0
#nu_1=nu_0
#rho_0=rho_1
#nu_0=nu_1
sigma_01=0.0#72.8e-3
#gravity
g=[0.0,0.0,-9.8]
useStokes=False
#
#domain
#
nd = 3
height = 1.0
length = 3.22
width  = 1.0
box_height = 0.161
box_width  = 0.403
box_length = 0.161
box_xy = [2.3955,.2985]
#
he = 0.5*box_width#2G
#he = 0.25*0.4*box_width#64*2G
#he = 0.5*0.25*0.4*box_width#64*2G

from boxInTank3dDomain import *
domain = boxInTank3d(L=[length,width,height],
                       box_xy=box_xy,
                       box_L=[box_length,box_width,box_height])
domain.writePoly("boxInFlume3d")
domain.writePLY("boxInFlume3d")
domain.writeAsymptote("boxInFlume3d")
#veg
from boxesInTank3dDomain import *
domain = boxesInTank3d(L=[length,width,height],
                       nBoxes_xy = [int(length/width)*4,4],
                       boxScale=[0.4,0.5,0.45])
domain.writePoly("boxesInFlume3d")
domain.writePLY("boxesInFlume3d")
domain.writeAsymptote("boxesInFlume3d")
bcsTimeDependent=False
#
#cylinder
#
# he = 0.2*box_width#2G
# #he = 0.25*0.4*box_width#64*2G
# #he = 0.5*0.25*0.4*box_width#64*2G
# from cylinderInTank3dDomain import *
# domain = cylinderInTank3d("cylinderInTank3d",
#                           L=[length,width,height],
#                           radius=0.1*width,
#                           center=(0.5*length,0.5*width,0.5*height),
#                           n_points_on_obstacle=2*5-2)
# domain.writePoly("cylinderInTank3d")
# domain.writePLY("cylinderInTank3d")
# domain.writeAsymptote("cylinderInTank3d")

boundaryTags = domain.boundaryTags

waterLevel = 0.5
Um=2.0

def shockSignedDistance(x):
    return x[2] - waterLevel

#
#time interval etc.
#
dt_init=he*0.001
T = 3.0
nDTout=300#500
runCFL = 0.33
#
#numerics
#
nLevels = 1
triangleOptions="VApq1.25q12ena%21.16e" % ((he**3)/6.0,)
print triangleOptions
applyCorrection=True
applyRedistancing=True
rdtimeIntegration='osher'
#rdtimeIntegration='newton'
freezeLevelSet=True
obstacleInFlume_quad_order = 3
useBackwardEuler=True
useBackwardEuler_ls=True
#subgrid error
lag_ns_subgridError=True
lag_ls_shockCapturing=True
#shock capturing diffusion
ns_shockCapturingFactor=0.1
ls_shockCapturingFactor=0.1
vof_shockCapturingFactor=0.1
rd_shockCapturingFactor=0.1
#epsilons for Heaviside/Dirac/etc smoothing
epsFact_density = 1.5
epsFact_viscosity = 1.5
epsFact_redistance = 0.33
epsFact_curvature=1.5
epsFact_consrv_heaviside=1.5
epsFact_consrv_dirac=1.5
epsFact_consrv_diffusion=10.0
epsFact_vof=1.5
#
usePETSc=True
spaceOrder=1
restrictFineSolutionToAllMeshes=False
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 1
