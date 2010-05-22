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
nu_1= 1.500e-5 * 1000.0#cek hack turbulence
sigma_01=0.0#72.8e-3
#gravity
g=[0.0,0.0,-9.8]
useStokes=False
#
#domain
#
nd = 3
#
genMesh=True
from vesselPoly3dDomain import *
domain = vesselPoly3d(fileprefix="Container06_test",vesselprefix="Container06_3d",outerBoxFactor=[10.0,10.0,0.1],offsetFactor=[0.6,0.0,0.0])
domain.writePoly("Container06_test")
length_scale = max(domain.L)

hull_beam = domain.hull_beam
hull_draft = domain.hull_draft
hull_length=  domain.hull_length
inflow_height= domain.L[2]
bottom_width = domain.L[1]
bottom_length= domain.L[0]
boundaryTags = domain.boundaryTags
print boundaryTags

height = domain.L[2] + domain.x[2]
waterLevel = 0.5*domain.L[2] + domain.x[2]#height

movingDomain=False#True
#
#residence time based on mean velocity
#
Um = 5.0#m/s 10 knots
RE = hull_length*Um/nu_0
Fr = Um/math.sqrt(math.fabs(g[2])*hull_length)
Frh = Um/math.sqrt(math.fabs(g[2])*waterLevel)
print "========================================REYNOLDS NUMBER = "+`RE`
print "========================================FROUDE(HULL LENGTH) NUMBER = "+`Fr`
print "========================================FROUDE(DEPTH) NUMBER = "+`Frh`
print "========================================SPEED[M/S] = %"+`Um`
residence_time = bottom_length/Um
#
#time interval etc.
#
T = 5.0#residence_time
print "========================================residence time = "+`residence_time`
print "========================================T = "+`T`
nDTout=50
runCFL = 1.5#0.33
#
#numerics
#
nLevels = 1
he = 0.25*hull_beam
dt_init= 0.001*he/Um
#triangleOptions="VApq1.25q13ena%f" % ((he**3)/6.0,)
#triangleOptions="KVApfenq1.25q12"#a%21.16e
triangleOptions="KVApfenq1.25q12"#a%21.16e
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
lag_ls_shockCapturing=True
#shock capturing diffusion
ns_shockCapturingFactor=0.9
ls_shockCapturingFactor=0.9
vof_shockCapturingFactor=0.9
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
