from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *
from proteus.Profiling import logEvent

from proteus import Context

#===============================================================================
# Context
#===============================================================================
ct = Context.Options([
    ("T", 3.0, "Time interval [0, T]"),
    ("onlySaveFinalSolution",False,"Only save the final solution"),
    ("Refinement",2, "Specify initial mesh size by giving number of cells in each direction"),
    ("spaceOrder",1,"FE space for velocity"),
    ("parallel",False,"Use parallel or not"),
    ("dt_fixed",0.005,"fixed time step"),
    ("nDTout",100,"output number"),
    ("use_sbm",True,"use sbm instead of imb"),
], mutable=True)


#===============================================================================
# Global parameters
#===============================================================================

if ct.use_sbm:
    USE_SBM=1
else:
    USE_SBM=0

Refinement = ct.Refinement
sedimentDynamics=False
genMesh = True
movingDomain = False
applyRedistancing = True
useOldPETSc = False
useSuperlu = True


parallel = ct.parallel
if parallel:
    usePETSc = True
    useSuperlu=False
else:
    usePETSc = False

timeDiscretization = 'vbdf'#vbdf'#'vbdf'  # 'vbdf', 'be', 'flcbdf'
spaceOrder = ct.spaceOrder
pspaceOrder = 1
useHex = False
useRBLES = 0.0
useMetrics = 1.0
applyCorrection = False
useVF = 0.0
useOnlyVF = False
useRANS = 0  # 0 -- None
             # 1 -- K-Epsilon
             # 2 -- K-Omega
openTop=True
fl_H = 0.41
# Input checks
if spaceOrder not in [1, 2]:
    print("INVALID: spaceOrder" + spaceOrder)
    sys.exit()

if useRBLES not in [0.0, 1.0]:
    print("INVALID: useRBLES" + useRBLES)
    sys.exit()

if useMetrics not in [0.0, 1.0]:
    print("INVALID: useMetrics")
    sys.exit()

#===============================================================================
# FE spaces
#===============================================================================
nd = 3

if spaceOrder == 1:
    hFactor = 1.0
    if useHex:
        basis = C0_AffineLinearOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd, 2)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd - 1, 2)
    else:
        basis = C0_AffineLinearOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd, 3)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd - 1, 3)
elif spaceOrder == 2:
    hFactor = 0.5
    if useHex:
        basis = C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd, 4)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd - 1, 4)
    else:
        basis = C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd, 5)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd - 1, 5)

if pspaceOrder == 1:
    if useHex:
        pbasis = C0_AffineLinearOnCubeWithNodalBasis
    else:
        pbasis = C0_AffineLinearOnSimplexWithNodalBasis
elif pspaceOrder == 2:
    if useHex:
        pbasis = C0_AffineLagrangeOnCubeWithNodalBasis
    else:
        pbasis = C0_AffineQuadraticOnSimplexWithNodalBasis

weak_bc_penalty_constant = 100.0
nLevels = 1
#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
structured = False

#===============================================================================
# Gauges
#===============================================================================

# new style PointGauges
# pointGauges = PointGauges(gauges = ((('u', 'v'), ((0.5, 0.5, 0), (1, 0.5, 0))), (('p',), ((0.5, 0.5, 0),))),
#                           activeTime=(0, 0.5),
#                           sampleRate=0,
#                           fileName='combined_gauge_0_0.5_sample_all.csv')

# lineGauges = LineGauges(gaugeEndpoints={'lineGauge_xtoH=0.825': ((0.495, 0.0, 0.0), (0.495, 1.8, 0.0))}, linePoints=20)
# #'lineGauge_x/H=1.653':((0.99,0.0,0.0),(0.99,1.8,0.0))
# lineGauges_phi = LineGauges_phi(lineGauges.endpoints, linePoints=20)

#===============================================================================
# Domain and mesh
#===============================================================================
L = [1.5, 0.41, 0.01]
he = old_div(L[0],(2**Refinement))
he*=0.5
he*=0.5

he = 0.01
L[2] = 4*he

import my_domain
# domain,boundaryTags = get_domain()
# domain,boundaryTags = my_domain.get_pseudo_3D_cylinder_domain(L=L,center=(0.2,0.2),radius=0.06,he=he)
domain,boundaryTags = my_domain.get_pseudo_3D_cylinder_box_domain(L=L,
                                                                 center=(0.2,0.2),
                                                                 radius=0.05,
                                                                 he=he,he2=0.5*he,
                                                                 x2=(0.2-0.06,0.2-0.06),x3=(0.2+0.3,0.2+0.06),
                                                                 )

domain.writePoly("mesh")
domain.writePLY("mesh")
domain.writeAsymptote("mesh")
#triangleOptions = "VApq30Dena%8.8f" % ((he ** 2) / 2.0,)
#triangleOptions = "KVApq10Dena"
triangleOptions="KVApq1.35q12feena"

logEvent("""Mesh generated using: tetgen -%s %s""" % (triangleOptions, domain.polyfile + ".poly"))

# boundaries = ['left', 'right', 'bottom', 'top', 'front', 'back']
# boundaryTags = dict([(key, i + 1) for (i, key) in enumerate(boundaries)])
# L = [1,1,1]
# domain = Domain.RectangularDomain(L)
# nnx=nny=nnz=3
#===============================================================================
# Time stepping
#===============================================================================
T=ct.T
dt_fixed = ct.dt_fixed#0.03
dt_init = 0.5*dt_fixed
runCFL=0.33
nDTout = ct.nDTout
# nDTout = int(T/dt_fixed)
dt_output = old_div(T,nDTout)
tnList = [0.0,dt_init]+[i*dt_output for i in range(1,nDTout+1)]

if ct.onlySaveFinalSolution == True:
    tnList = [0.0,dt_init,ct.T]

#===============================================================================
# Numerical parameters
#===============================================================================
ns_forceStrongDirichlet = True
ns_sed_forceStrongDirichlet = False
if useMetrics:
    ns_shockCapturingFactor  = 0.5
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.5
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.5
    vof_shockCapturingFactor = 0.5
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.5
    rd_shockCapturingFactor  = 0.5
    rd_lag_shockCapturing = False
    epsFact_density    = 1.5
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = True
    kappa_shockCapturingFactor = 0.25
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.25
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0
else:
    ns_shockCapturingFactor = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ns_sed_shockCapturingFactor = 0.9
    ns_sed_lag_shockCapturing = True
    ns_sed_lag_subgridError = True
    ls_shockCapturingFactor = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref = 1.0
    ls_sc_beta = 1.0
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    vos_shockCapturingFactor = 0.9
    vos_lag_shockCapturing = True
    vos_sc_uref = 1.0
    vos_sc_beta = 1.0
    rd_shockCapturingFactor = 0.9
    rd_lag_shockCapturing = False
    epsFact_density = 1.5
    epsFact_viscosity = epsFact_curvature = epsFact_vof = epsFact_vos = epsFact_consrv_heaviside = epsFact_consrv_dirac = \
        epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 0.1
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True  #False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True  #False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0

ns_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
ns_sed_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
vof_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
vos_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
ls_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
rd_nl_atol_res = max(1.0e-10, 0.05 * he)
mcorr_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
kappa_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
dissipation_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
phi_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
pressure_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)

#===============================================================================
# turbulence
#===============================================================================
ns_closure = 0  #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
ns_sed_closure = 0  #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4
#===============================================================================
# Physics parameters
#===============================================================================
# Water
rho_0 = 1.0e0
nu_0 = 1.0e-3

# Air
rho_1 = rho_0#1.205
nu_1 = nu_0#1.500e-5

# Sediment

rho_s = rho_0
nu_s = 10000.0*nu_0
dragAlpha = 0.0

# Surface tension
sigma_01 = 0.0

# Gravity
g = [0.0, 0.0, 0.0]

# Initial condition
waterLine_x = 0.75
waterLine_z = 1.6

#===============================================================================
# Global functions
#===============================================================================
U = 1.5 # this is the inlet max velocity not the mean velocity
def velRamp(t,x):
#     return U*16.0*x[1]*(fl_H-x[1])/fl_H**2*x[2]*(L[2]-x[2])/L[2]**2
    return U*4*x[1]*(fl_H-x[1])/fl_H**2

#===============================================================================
# Particles
#===============================================================================
def particle_sdf(t, x):
    cx = 0.2
    cy = 0.2
    r = math.sqrt( (x[0]-cx)**2 + (x[1]-cy)**2)
    n = (old_div((x[0]-cx),(r+1e-10)),old_div((x[1]-cy),(r+1e-10)), 0.0)
    return  r - 0.05,n

def particle_vel(t, x):
    return (0.0,0.0,0.0)
