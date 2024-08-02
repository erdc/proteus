from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *
from proteus.Profiling import logEvent
from .parameters import *

manufactured_solution = 2 #1: u.n!=0, 2: u.n=0

# ----- PARAMETERS ABOUT STABILIZATION OF NS ----- #
# see parameters.py

# ----- OTHER IMPORTANT PARAMETERS FOR CONVERGENCE TEST ----- #
KILL_PRESSURE_TERM = False #True just for debigging
Refinement = 0
mu = 1.0 
ns_forceStrongDirichlet = False

# ----- ENTROPY VISCOSITY PARAMETERS ----- #
cE = 1.0
cMax = 1.0

# ----- numerical options to RANS3PF ----- #
fixNullSpace_PresInc = True
INTEGRATE_BY_PARTS_DIV_U_PresInc = True
CORRECT_VELOCITY = True

# ----- Discretization -- input options ----- #
genMesh = True
movingDomain = False
useOldPETSc = False
useSuperlu = True
timeDiscretization = 'vbdf' # 'vbdf', 'be', 'flcbdf'
spaceOrder = 2
pspaceOrder = 1
useHex = False
useMetrics = 0.0
useVF = 0.0
useRBLES=0.0

# Input checks
if spaceOrder not in [1, 2]:
    print("INVALID: spaceOrder" + spaceOrder)
    sys.exit()

if useMetrics not in [0.0, 1.0]:
    print("INVALID: useMetrics")
    sys.exit()

#  Discretization
nd = 2

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

# Domain and mesh
L = (1.0, 1.0)
he = L[0]/float(4*Refinement-1)
he*=0.5
he*=0.5

weak_bc_penalty_constant = 1E6
nLevels = 1
#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
structured = False

if useHex:
    nnx = 4 * Refinement + 1
    nny = 2 * Refinement + 1
    hex = True
    domain = Domain.RectangularDomain(L)
else:
    boundaries = ['left', 'right', 'bottom', 'top', 'front', 'back']
    boundaryTags = dict([(key, i + 1) for (i, key) in enumerate(boundaries)])
    if structured:
        nnx = 4 * Refinement
        nny = 4 * Refinement
    else:
        vertices = [[0.0, 0.0],  #0
                    [L[0], 0.0],  #1
                    [L[0], L[1]],  #2
                    [0.0, L[1]]]  #3
        vertexFlags = [boundaryTags['bottom'],
                       boundaryTags['bottom'],
                       boundaryTags['top'],
                       boundaryTags['top']]
        segments = [[0, 1],
                    [1, 2],
                    [2, 3],
                    [3, 0]]
        segmentFlags = [boundaryTags['bottom'],
                        boundaryTags['right'],
                        boundaryTags['top'],
                        boundaryTags['left']]
        regions = [[1.2, 0.6]]
        regionFlags = [1]

        domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                      vertexFlags=vertexFlags,
                                                      segments=segments,
                                                      segmentFlags=segmentFlags,
                                                      regions=regions,
                                                      regionFlags=regionFlags)
        #go ahead and add a boundary tags member
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        triangleOptions = "VApq30Dena%8.8f" % ((he**2)/2.0,)

domain.MeshOptions.nn = nn
domain.MeshOptions.nnx = nnx
domain.MeshOptions.nny = nny
domain.MeshOptions.nnz = nnz
domain.MeshOptions.triangleFlag=0
domain.MeshOptions.triangleOptions=triangleOptions

# Time stepping
T=1.0
dt_fixed = 0.1
dt_init = min(0.1*dt_fixed,0.001)
runCFL=0.33
nDTout = int(round(T/dt_fixed))

# Numerical parameters
ns_sed_forceStrongDirichlet = False
if useMetrics:
    ns_shockCapturingFactor  = 0.5
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.5
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    vof_shockCapturingFactor = 0.5
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    rd_shockCapturingFactor  = 0.5
    rd_lag_shockCapturing = False
    epsFact_density    = 3.0
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

#turbulence
ns_closure = 0  #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
ns_sed_closure = 0  #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
    
# Water. Fake parameters for convergence test
rho_0 = 1.0
nu_0 = 1.0

# Air
rho_1 = 1.0
nu_1 = 1.0

# Sediment
rho_s = rho_0
nu_s = 10000.0*nu_0
dragAlpha = 0.0

# Surface tension
sigma_01 = 0.0

# Gravity
g = [0.0, 0.0]

