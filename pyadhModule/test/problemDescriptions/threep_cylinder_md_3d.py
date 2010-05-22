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
nu_1= 1.500e-5
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
inflow_height=0.41
he = 0.1*inflow_height
bottom_width=he
bottom_length=inflow_height
cylinder_radius=inflow_height/10.0
cylinder_center = (bottom_length/2.0,bottom_width/2.0,inflow_height/2.0)

from cylinder3dDomain import *
domain = cylinder3d(fileprefix="cylinder3d",
                    cross_section=circular_cross_section,
                    height=inflow_height,
                    length=bottom_length,
                    width=bottom_width,
                    radius=cylinder_radius,
                    center = cylinder_center,
                    n_points_on_obstacle=2*5-2)
domain.writePoly("cylinder3d")
domain.writeAsymptote("cylinder3d")
boundaryTags = domain.boundaryTags
waterLevel = 0.5*inflow_height
useShock=False#True
movingDomain=False#True
#
#residence time based on mean velocity
#
RE = 1.0
Um = nu_0*RE/(2.0*cylinder_radius)
#RE = 2.0*cylinder_radius*Um/nu_0
Profiling.logEvent("REYNOLDS NUMBER = "+`RE`)
residence_time = bottom_length/Um
#Um=0.0
#
#time interval etc.
#
dt_init=0.01*residence_time
T = 2#10.0*residence_time
nDTout=100
runCFL = 1.5
#
#numerics
#
nLevels = 1
triangleOptions="VpAq1.25ena%f" % ((he**3)/6.0,)

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
rd_shockCapturingFactor=0.9
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
parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.element
nLayersOfOverlapForParallel = 1
