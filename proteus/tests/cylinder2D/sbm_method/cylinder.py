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

ct = Context.Options([
    ("T", 4.0, "Time interval [0, T]"),
    ("he",0.04, "maximum size of edges"),
    ("onlySaveFinalSolution",False,"Only save the final solution"),
    ("vspaceOrder",2,"FE space for velocity"),
    ("pspaceOrder",1,"FE space for pressure"),
    ("use_sbm",True,"use sbm instead of imb"),
    ("use_regions",False, "use refinement regions")
], mutable=True)


if ct.use_sbm:
    USE_SBM=1
else:
    USE_SBM=0

#  Discretization -- input options
sedimentDynamics=False
genMesh = True
movingDomain = False
applyRedistancing = True
useOldPETSc = False
useSuperlu = True
timeDiscretization = 'vbdf'#vbdf'#'vbdf'  # 'vbdf', 'be', 'flcbdf'
spaceOrder = ct.vspaceOrder
pspaceOrder = ct.pspaceOrder
useHex = False
useRBLES = 0.0
useMetrics = 1.0
applyCorrection = True
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
#L = (0.584,0.350)
L = (2.2, 0.41)
he = ct.he
weak_bc_penalty_constant = 100.0
nLevels = 1
#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
structured = False



# new style PointGauges
# pointGauges = PointGauges(gauges = ((('u', 'v'), ((0.5, 0.5, 0), (1, 0.5, 0))), (('p',), ((0.5, 0.5, 0),))),
#                           activeTime=(0, 0.5),
#                           sampleRate=0,
#                           fileName='combined_gauge_0_0.5_sample_all.csv')

# lineGauges = LineGauges(gaugeEndpoints={'lineGauge_xtoH=0.825': ((0.495, 0.0, 0.0), (0.495, 1.8, 0.0))}, linePoints=20)
# #'lineGauge_x/H=1.653':((0.99,0.0,0.0),(0.99,1.8,0.0))
# lineGauges_phi = LineGauges_phi(lineGauges.endpoints, linePoints=20)

if useHex:
    nnx = ceil(L[0]/ct.he)
    nny = ceil(L[1]/ct.he)
    hex = True
    domain = Domain.RectangularDomain(L)
else:
    boundaries = ['bottom', 'right', 'top', 'left', 'front', 'back']
    boundaryTags = dict([(key, i + 1) for (i, key) in enumerate(boundaries)])
    if structured:
        nnx = ceil(L[0]/ct.he)
        nny = ceil(L[1]/ct.he)
    else:
        if ct.use_regions:
            vertices = [[0.0, 0.0],  #0
                        [L[0], 0.0],  #1
                        [L[0], L[1]],  #2
                        [0.0, L[1]],  #3
                        [0.2-0.16,L[1]*0.2],
                        [0.2-0.16,L[1]*0.8],
                        [0.2+0.3,L[1]*0.8],
                        [0.2+0.3,L[1]*0.2],
                        # the following are set for refining the mesh
                    [0.2-0.06,0.2-0.06],
                        [0.2-0.06,0.2+0.06],
                        [0.2+0.06,0.2+0.06],
                        [0.2+0.06,0.2-0.06]]

                    
                    
            vertexFlags = [boundaryTags['bottom'],
                           boundaryTags['bottom'],
                           boundaryTags['top'],
                           boundaryTags['top'],
                           # the interior vertices should be flaged to 0
                           0, 0, 0, 0,
                           0, 0, 0, 0 ]
            
            segments = [[0, 1],
                        [1, 2],
                        [2, 3],
                        [3, 0],
                        #Interior segments
                        [4, 5],
                        [5, 6],
                        [6, 7],
                        [7,4],
                        [8,9],
                        [9,10],
                        [10,11],
                        [11,8]]
            segmentFlags = [boundaryTags['bottom'],
                            boundaryTags['right'],
                            boundaryTags['top'],
                            boundaryTags['left'],
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0]
            
            regions = [[0.95*L[0], 0.2],[0.2-0.15,0.2],[0.2,0.2]]
            regionFlags = [1,2,3]
            regionConstraints=[0.5*he**2,0.5*(old_div(he,2.0))**2,0.5*(old_div(he,6.0))**2]
            #        for gaugeName,gaugeCoordinates in pointGauges.locations.iteritems():
            #            vertices.append(gaugeCoordinates)
            #            vertexFlags.append(pointGauges.flags[gaugeName])
            
            # for gaugeName, gaugeLines in lineGauges.linepoints.iteritems():
            #     for gaugeCoordinates in gaugeLines:
            #         vertices.append(gaugeCoordinates)
            #         vertexFlags.append(lineGauges.flags[gaugeName])
            domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                          vertexFlags=vertexFlags,
                                                          segments=segments,
                                                          segmentFlags=segmentFlags,
                                                          regions=regions,
                                                          regionFlags=regionFlags,
                                                          regionConstraints=regionConstraints)
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
            
            regions = [[0.95*L[0], 0.2]]
            regionFlags = [1]
            regionConstraints=[0.5*he**2]
            #        for gaugeName,gaugeCoordinates in pointGauges.locations.iteritems():
            #            vertices.append(gaugeCoordinates)
            #            vertexFlags.append(pointGauges.flags[gaugeName])
            
            # for gaugeName, gaugeLines in lineGauges.linepoints.iteritems():
            #     for gaugeCoordinates in gaugeLines:
            #         vertices.append(gaugeCoordinates)
            #         vertexFlags.append(lineGauges.flags[gaugeName])
            domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                          vertexFlags=vertexFlags,
                                                          segments=segments,
                                                          segmentFlags=segmentFlags,
                                                          regions=regions,
                                                          regionFlags=regionFlags,
                                                          regionConstraints=regionConstraints)
            
        #go ahead and add a boundary tags member
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        #triangleOptions = "VApq30ena%8.8f" % ((he ** 2) / 2.0,)
        triangleOptions = "VApq30Dena"

logEvent("""Mesh generated using: tetgen -%s %s""" % (triangleOptions, domain.polyfile + ".poly"))
# Time stepping
T=ct.T
dt_fixed = 0.005#0.03
dt_init = 0.005#min(0.1*dt_fixed,0.001)
runCFL=0.33
nDTout = int(round(old_div(T,dt_fixed)))
tnList = [0.0,dt_init]+[i*dt_fixed for i in range(1,nDTout+1)]

if ct.onlySaveFinalSolution == True:
    tnList = [0.0,dt_init,ct.T]
# Numerical parameters
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

#turbulence
ns_closure = 0  #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
ns_sed_closure = 0  #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4
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
g = [0.0, 0.0]

# Initial condition
waterLine_x = 0.75
waterLine_z = 1.6

U = 1.5 # this is the inlet max velocity not the mean velocity
def velRamp(t):
    if t < 2.0:
        return t*U/2.0
    else:
        return U
    return U



def signedDistance(x):
    return x[1]-old_div(L[1],2)

def particle_sdf(t, x):
    cx = 0.2
    cy = 0.2
    r = math.sqrt( (x[0]-cx)**2 + (x[1]-cy)**2)
    n = (old_div((x[0]-cx),r),old_div((x[1]-cy),r))
    return  r - 0.05,n

def particle_vel(t, x):
    return (0.0,0.0)

