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
    ("he",0.05,"max element diameter"),
    ("schur_solver", 'two_phase_PCD', "preconditioner type"),
    ("ns_forceStrongDirichlet",True,"boundary condition type"),
    ("boundary_condition_type",'fs',"Free or no slip boundary"),
    ("reynolds_number",10,"Simulation's Reynolds number"),
    ])

nd = 2
spaceOrder=1
useVF = 0.0
usePETSc = True
useOnlyVF = True
useMetrics = True
he = opts.he

# Input checks
if spaceOrder not in [1,2]:
    print("INVALID: spaceOrder" + spaceOrder)
    sys.exit()
if spaceOrder == 1:
    hFactor=1.0
    basis=C0_AffineLinearOnSimplexWithNodalBasis
    elementQuadrature = SimplexGaussQuadrature(nd,5)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,5)
elif spaceOrder == 2:
    hFactor=0.5
    basis=C0_AffineQuadraticOnSimplexWithNodalBasis
    elementQuadrature = SimplexGaussQuadrature(nd,5)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,5)

nLevels = 1
L=(6.0,1.0)
x0=(-1.0,0.0)
boundaries = ['left', 'right', 'bottom', 'top', 'front', 'back']
boundaryTags = dict([(key, i + 1) for (i, key) in enumerate(boundaries)])
vertices    = [[-1.0,0.0],[-1.0,1.0],[5.0,1.0],[5.0,-1.0],[0.0,-1.0],[0.0,0.0]]
vertexFlags = [boundaryTags['left'],
               boundaryTags['left'],
               boundaryTags['right'],
               boundaryTags['bottom'],
               boundaryTags['bottom'],
               boundaryTags['bottom']]
segments = [[0, 1], [1, 2], [2, 3], [3, 4],[4,5], [5,0]]
segmentFlags = [boundaryTags['left'],
                boundaryTags['top'],
                boundaryTags['right'],
                boundaryTags['bottom'],
                boundaryTags['bottom'],
                boundaryTags['bottom']]
regions = [[1.0, 1.0]]
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
domain.MeshOptions.triangleOptions = "VApq30Dena%8.8f" % ((he ** 2)/2.0,)

logEvent("""Mesh generated using: tetgen -%s %s""" % (triangleOptions, domain.polyfile + ".poly"))




# Numerical parameters

ns_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
vof_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
ls_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
rd_nl_atol_res = max(1.0e-10, 0.005 * he)
mcorr_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
kappa_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)
dissipation_nl_atol_res = max(1.0e-10, 0.001 * he ** 2)

#Single-phase example with Re=10
# rho_0 = 1.0
# nu_0 = 1/5.
# rho_1 = 1.0
# nu_1 = 1/5.

#Single-phase example with Re=100
# rho_0 = 1.0
# nu_0 = 1.0/ 50.
# rho_1 = 1.0
# nu_1 = 1.0 / 50.

#Two-phase example with Re=10
# rho_0 = 1.0
# nu_0 = 1.0 / 5.
# rho_1 = 1.205/998.2
# nu_1 = (1.0 / 5.) * ( 1.500e-5 / 1.004e-6)

#Two-phase example with Re=100
assert opts.reynolds_number in [100,10]
if opts.reynolds_number == 100:
    rho_0 = 1.0 #1.205 / 998.2
    nu_0 = 1.0/50.
    rho_1 =  1.205/998.2
    nu_1 = (1.0/50.) * ( 1.500e-5/1.004e-6)

elif opts.reynolds_number == 10:
    rho_0 = 1.0
    nu_0 = 1.0/5.
    rho_1 = 1.205/998.2
    nu_1 = (1.0/5.) * (1.500e-5/1.004e-6)
# Gravity
g = [0.0,0.0]


ns_forceStrongDirichlet = opts.ns_forceStrongDirichlet
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
epsFact_density = 0.0
epsFact_viscosity = epsFact_curvature = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = \
    epsFact_density
epsFact_redistance = 0.33
epsFact_consrv_diffusion = 0.1
redist_Newton = True
kappa_shockCapturingFactor = 0.25
kappa_lag_shockCapturing = True
kappa_sc_uref = 1.0
kappa_sc_beta = 1.0
dissipation_shockCapturingFactor = 0.25
dissipation_lag_shockCapturing = True
dissipation_sc_uref = 1.0
dissipation_sc_beta = 1.0
