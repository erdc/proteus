from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *
from proteus.Profiling import logEvent

#  Discretization -- input options
#Refinement = 20#45min on a single core for spaceOrder=1, useHex=False
Refinement = 5 #45min on a single core for spaceOrder=1, useHex=False
genMesh=True
useOldPETSc=False
useSuperlu=False
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
applyCorrection=True
useVF = 1.0
useOnlyVF = True
redist_Newton = True
useRANS = 0 # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega
# Input checks
if spaceOrder not in [1,2]:
    print "INVALID: spaceOrder" + spaceOrder
    sys.exit()

if useRBLES not in [0.0, 1.0]:
    print "INVALID: useRBLES" + useRBLES
    sys.exit()

if useMetrics not in [0.0, 1.0]:
    print "INVALID: useMetrics"
    sys.exit()

#  Discretization
nd = 3
if spaceOrder == 1:
    hFactor=1.0
    if useHex:
	 basis=C0_AffineLinearOnCubeWithNodalBasis
         elementQuadrature = CubeGaussQuadrature(nd,2)
         elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,2)
    else:
    	 basis=C0_AffineLinearOnSimplexWithNodalBasis
         elementQuadrature = SimplexGaussQuadrature(nd,3)
         elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
elif spaceOrder == 2:
    hFactor=0.5
    if useHex:
	basis=C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd,4)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,4)
    else:
	basis=C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd,4)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

# Domain and mesh
L = (0.8,0.2,0.2)
#he = L[0]/float(4*Refinement-1)
#he*=0.5
he = L[0]/12.0
he*=0.5#128
#he*=0.5#1024
quasi2D = False#True

if quasi2D:
    L = (L[0],he,L[2])

# Domain and mesh
nLevels = 1
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
structured=False
if useHex:
    nnx=4*Refinement+1
    nny=1*Refinement+1
    nnz=2*Refinement+1
    if quas2D:
        nny=2
    hex=True
    domain = Domain.RectangularDomain(L)
else:
    boundaries=['left','right','bottom','top','front','back']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    if structured:
        nnx=4*Refinement
        nny=2*Refinement
    else:
        vertices=[[0.0,0.0,0.0],#0
                  [L[0],0.0,0.0],#1
                  [L[0],L[1],0.0],#2
                  [0.0,L[1],0.0],#3
                  [0.0,0.0,L[2]],#4
                  [L[0],0.0,L[2]],#5
                  [L[0],L[1],L[2]],#6
                  [0.0,L[1],L[2]]]#7
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
        regions=[[0.5*L[0],0.5*L[1],0.5*L[2]]]
        regionFlags=[1]
        domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                     vertexFlags=vertexFlags,
                                                     facets=facets,
                                                     facetFlags=facetFlags,
                                                     regions=regions,
                                                     regionFlags=regionFlags)
        #go ahead and add a boundary tags member
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        domain.numAdaptSteps=1
        triangleOptions="VApq1.4q12feena%21.16e" % ((he**3)/6.0,)

logEvent("""Mesh generated using: tetgen -%s %s"""  % (triangleOptions,domain.polyfile+".poly"))
# Time stepping
T=0.03
dt_fixed = 0.01
dt_init = min(0.1*dt_fixed,0.001)
runCFL=0.33
nDTout = int(round(T/dt_fixed))

# Numerical parameters
ns_forceStrongDirichlet = True
if useMetrics:
    ns_shockCapturingFactor  = 0.9
    ns_lag_shockCapturing = False
    ns_lag_subgridError = False
    ls_shockCapturingFactor  = 0.9
    ls_lag_shockCapturing = False
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.5
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = False
    vof_sc_uref = 1.0
    vof_sc_beta = 1.5
    rd_shockCapturingFactor  = 0.9
    rd_lag_shockCapturing = False
    epsFact_density    = 1.5
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False
else:
    ns_shockCapturingFactor  = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref  = 1.0
    vof_sc_beta  = 1.0
    rd_shockCapturingFactor  = 0.9
    rd_lag_shockCapturing = False
    epsFact_density    = 1.5
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

ns_nl_atol_res = max(1.0e-6,0.000001*he**2)
vof_nl_atol_res = max(1.0e-6,0.000001*he**2)
ls_nl_atol_res = max(1.0e-6,0.000001*he**2)
rd_nl_atol_res = max(1.0e-6,0.000001*he**2)
mcorr_nl_atol_res = max(1.0e-6,0.000001*he**2)
kappa_nl_atol_res = max(1.0e-6,0.000001*he**2)
dissipation_nl_atol_res = max(1.0e-6,0.000001*he**2)

#turbulence
ns_closure=2 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4
# Water
rho_0 = 1000.0
nu_0  = 1.0e-6

# Air
rho_1 = 1000.0
nu_1  = 1.0e-6

# Surface tension
sigma_01 = 0.0

# Gravity
g = [0.0,0.0,0.0]

# Initial condition
waterLine_x = 0.40
waterLine_z = 0.0

def signedDistance(x):
	return (x[0]-waterLine_x)
#    phi_x = x[0]-waterLine_x
#    phi_z = x[2]-waterLine_z
#    if phi_x < 0.0:
#        if phi_z < 0.0:
#            return max(phi_x,phi_z)
#        else:
#            return phi_z
#    else:
#        if phi_z < 0.0:
#            return phi_x
#        else:
#            return sqrt(phi_x**2 + phi_z**2)
