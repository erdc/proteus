"""
A helper module for doing air/water flow around a moving rigid cylinder in 2D
"""
from math import *
from pyadh import *

#
#material properties
#
#water
rho_0=998.2
nu_0=1.004e-6
#air
rho_1=1.205
nu_1= 1.500e-5# *1.0e5#cek hack damp turbulence in air phase
sigma_01=0.0#72.8e-3
#gravity
g=[0.0,0.0,-9.8]
useStokes=False
#
#domain
#
nd = 3
hull_beam = 0.238
hull_draft = 0.095
hull_length=  1.905
height=2.0*hull_draft#2.0*(3.0*hull_draft)
length=2.0*hull_length#
#length = 4.0*hull_length+hull_length
width=4.0*hull_length+hull_beam
width=3.0*hull_beam
height=3.0*hull_draft
length=2.0*hull_length
#width=4.0*hull_beam#4.0*hull_length+hull_beam
#
#width*=8.0
#length*=4.0
#
#he = 0.5*hull_draft#length
he = hull_draft#length
width=2.0*hull_beam
#he = 0.25*hull_draft#length (64 procs)
n_points_draft=max(2,int(ceil(hull_draft/he))) + 1
n_points_length=max(2,int(ceil(hull_length/he))) + 1

from wigley3dDomain import *
domain = wigley3d("wigley3d",
                  height,
                  length,
                  width,
                  hull_beam,
                  hull_draft,
                  hull_length,
                  (0.65*hull_length,0.5*width,0.5*height),
                  n_points_draft,
                  n_points_length)
domain.writePoly("wigley3d")
domain.writePLY("wigley3d")
domain.writeAsymptote("wigley3d")
boundaryTags = domain.boundaryTags
print boundaryTags
waterLevel = 0.5*height
openTop = True
openSides = False#True
smoothBottom = False
smoothObstacle = False
rampInitialConditions = False
#
# openTop = False
# openSides = False
# smoothBottom = True
# smoothObstacle = False
movingDomain=False#True
#
#residence time based on mean velocity
#
Um = 5.0#m/s 10 knots
RE = hull_beam*Um/nu_0
#RE = 4.2e6
#Um = RE*nu_0/hull_length
Fr = Um/math.sqrt(math.fabs(g[2])*hull_length)
Frh = Um/math.sqrt(math.fabs(g[2])*waterLevel)
print "========================================REYNOLDS NUMBER = "+`RE`
print "========================================FROUDE(HULL LENGTH) NUMBER = "+`Fr`
print "========================================FROUDE(DEPTH) NUMBER = "+`Frh`
print "========================================SPEED[M/S] = %"+`Um`
residence_time = length/Um
#
#time interval etc.
#
dt_init=1.0e-2#0.001*residence_time
T = 5.0*residence_time
print "========================================residence time = "+`residence_time`
print "========================================T = "+`T`
nDTout=100
runCFL = 0.33
#
#numerics
#
nLevels = 1
triangleOptions="VApq1.25q13fena%f" % ((he**3)/6.0,)
print triangleOptions

applyCorrection=True
applyRedistancing=True
rdtimeIntegration='osher'
freezeLevelSet=True
quad_order = 3
useBackwardEuler=True
useBackwardEuler_ls=True
#subgrid error
lag_ns_subgridError=True
lag_ns_shockCapturing=True
lag_ls_shockCapturing=True
#shock capturing diffusion
ns_shockCapturingFactor=0.33
ls_shockCapturingFactor=0.33
vof_shockCapturingFactor=0.33
rd_shockCapturingFactor=0.33
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
