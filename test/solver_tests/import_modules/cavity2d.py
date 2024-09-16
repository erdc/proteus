from math import *
import proteus.MeshTools
from proteus import Domain
from proteus.default_n import *
from proteus.default_so import *
from proteus import defaults
defaults.reset_default_n()
defaults.reset_default_so()
from proteus import Context

opts=Context.Options([
    ("he",0.2,"max element diameter"),
    ("schur_solver", 'two_phase_PCD', "preconditioner type"),
    ])

nd = 2
spaceOrder=1
Refinement=1
useHex=False
useVF = 0.0
#points_on_grain = 2
DX = opts.he #2 #0.02 #0.08, 0.04, 0.02
usePETSc = True
useOnlyVF = True
useMetrics = True
useOldPETSc = False
useSuperlu = False
timeDiscretization = 'be'

structured = False
L= (-1.0, 1.0)
#he = 0.1
he = opts.he

# Input checks
if spaceOrder not in [1,2]:
    print("INVALID: spaceOrder" + spaceOrder)
    sys.exit() 
if spaceOrder == 1:
    hFactor=1.0
    if useHex:
         basis=C0_AffineLinearOnCubeWithNodalBasis
         elementQuadrature = CubeGaussQuadrature(nd,5)
         elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,5)         
    else:
         basis=C0_AffineLinearOnSimplexWithNodalBasis
         elementQuadrature = SimplexGaussQuadrature(nd,5)
         elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,5)         
elif spaceOrder == 2:
    hFactor=0.5
    if useHex:    
        basis=C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd,5)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,5)    
    else:    
        basis=C0_AffineQuadraticOnSimplexWithNodalBasis 
        elementQuadrature = SimplexGaussQuadrature(nd,5)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,5)

nLevels = 1

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
        nny = 2 * Refinement
        domain = Domain.RectangularDomain(L)
    else:
        vertices = [[L[0], L[0]],  #0
                    [L[1], L[0]],  #1
                    [L[1], L[1]],  #2
                    [L[0], L[1]]]  #3
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
        # ARB - for two phase problem will need to change this...
        regions = [[1.0, 1.0]]
        regionFlags = [1]

        #        for gaugeName,gaugeCoordinates in pointGauges.locations.iteritems():
        #            vertices.append(gaugeCoordinates)
        #            vertexFlags.append(pointGauges.flags[gaugeName])

#        for gaugeName, gaugeLines in lineGauges.linepoints.iteritems():
#            for gaugeCoordinates in gaugeLines:
#                vertices.append(gaugeCoordinates)
#                vertexFlags.append(lineGauges.flags[gaugeName])
        domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                      vertexFlags=vertexFlags,
                                                      segments=segments,
                                                      segmentFlags=segmentFlags,
                                                      regions=regions,
                                                      regionFlags=regionFlags)
        #go ahead and add a boundary tags member
        domain.boundaryTags = boundaryTags
        #domain.writePoly("mesh")
        #domain.writePLY("mesh")
        #domain.writeAsymptote("mesh")
        domain.MeshOptions.triangleOptions = "VApq30Dena%8.8f" % ((he**2)/2.0,)

        #logEvent("""Mesh generated using: tetgen -%s %s""" % (triangleOptions, domain.polyfile + ".poly"))

domain.MeshOptions.nnx = nnx
domain.MeshOptions.nny = nny

#boundaryTags=domain.boundaryFlags
# Time stepping
#T = 1.0
# T=  8.00e0
# runCFL = 0.9
# dt_fixed = 1.0e-2  #10.0 #1.0e-2 #2.5e-1
# dt_init = 1.0e-3
# nDTout = int(T/dt_fixed)
# dt_init = min(dt_init,0.5*dt_fixed)

useBackwardEuler = False

# Numerical parameters

ns_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
vof_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
ls_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
rd_nl_atol_res = max(1.0e-10, 0.005 * he)
mcorr_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
kappa_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
dissipation_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)

# Reynolds number = ?
rho_0 = 1.205
nu_0 = 1.500e-5/1.500e-5
rho_1 = 998.2
nu_1 = 1.004e-6/1.500e-5

# Gravity
g = [0.0,0.0]

# ToDo (ARB) - I think this is needed for unstructured meshes so
# I've added a conditional statement.  You may need to confirm at some
# point.
if structured is False:
    triangleOptions= "pAq30.0Dena%f" % (.5*DX**2)  #% (0.5*(DX)**2,)
    print(triangleOptions)
    #genMesh=True
    #domain.writePLY('cavity2D')
    #domain.writePoly('cavity2D')
    domain.polyfile=os.path.dirname(os.path.abspath(__file__))+"/../"+"cavity2D"
    domain.MeshOptions.genMesh=False 
    domain.MeshOptions.triangleOptions = triangleOptions


ns_forceStrongDirichlet = True
if useMetrics:
    ns_shockCapturingFactor = 0.0
    ns_lag_shockCapturing = False
    ns_lag_subgridError = False
    ls_shockCapturingFactor = 0.5
    ls_lag_shockCapturing = True
    ls_sc_uref = 1.0
    ls_sc_beta = 1.0
    vof_shockCapturingFactor = 0.5
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    rd_shockCapturingFactor = 0.5
    rd_lag_shockCapturing = False
    epsFact_density = 1.5
    epsFact_viscosity = epsFact_curvature = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = \
        epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 0.1
    redist_Newton = True
    kappa_shockCapturingFactor = 0.25
    kappa_lag_shockCapturing = True  #False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.25
    dissipation_lag_shockCapturing = True  #False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0
else:
    ns_shockCapturingFactor = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref = 1.0
    ls_sc_beta = 1.0
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    rd_shockCapturingFactor = 0.9
    rd_lag_shockCapturing = False
    epsFact_density = 1.5
    epsFact_viscosity = epsFact_curvature = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = \
        epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 1.0
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True  #False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True  #False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0
