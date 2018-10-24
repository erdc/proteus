from __future__ import print_function
from __future__ import division
from past.utils import old_div
import numpy as np
from proteus import Domain
from proteus.mprans import SpatialTools as st
from proteus.mbd import ChRigidBody as crb

rho_0 = 1000.
nu_0 = 1.004e-6
rho_1 = 1.205
nu_1 = 1.500e-5
sigma_01 = 0.0
g = [0., -9.81]
he = 0.2
water_level = 2.5

# GEOMETRY

domain = Domain.PlanarStraightLineGraphDomain()

tank_dim = [5.,5.]
tank = st.Tank2D(domain, dim=tank_dim)
rect = st.Rectangle(domain, dim=[1.,1.], coords=[old_div(tank_dim[0],2.), old_div(tank_dim[1],2.)])
rect.setHoles(holes=np.array([rect.coords]))

domain.MeshOptions.he = he

# BOUNDARY CONDITIONS

tank.BC['x+'].setNoSlip()
tank.BC['x-'].setNoSlip()
tank.BC['y-'].setNoSlip()
tank.BC['y+'].setAtmosphere()

rect.BC['x+'].setNoSlip()
rect.BC['x-'].setNoSlip()
rect.BC['y+'].setNoSlip()
rect.BC['y-'].setNoSlip()

# CHRONO

system = crb.ProtChSystem(gravity=np.array([0.,-9.81,0.]))
body = crb.ProtChBody(system=system)
body.attachShape(rect)
body.ChBody.SetMass(500.)
body.ChBody.SetBodyFixed(True)  # fixing body

# OTHER PARAMS
st.assembleDomain(domain)

max_flag = max(domain.vertexFlags+domain.segmentFlags+domain.facetFlags)
flags_rigidbody = np.zeros(max_flag+1, dtype='int32')
for key in rect.boundaryTags_global:
    flags_rigidbody[rect.boundaryTags_global[key]] = 1


from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus.default_n import *

addedMass = True
movingDomain=True
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = old_div(10.,nu_0)
dt_init = 0.001
dt_fixed = 0.001
T = 0.002
nDTout = 1
timeIntegration = "backwardEuler"
runCFL = 0.4
nsave = 0

#----------------------------------------------------

#  Discretization -- input options
useOldPETSc=False
useSuperlu = not True
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
useVF = 1.0
useOnlyVF = False
useRANS = 0.
            # 1 -- K-Epsilon
            # 2 -- K-Omega, 1998
            # 3 -- K-Omega, 1988
# Input checks
if spaceOrder not in [1,2]:
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
    hFactor=1.0
    if useHex:
	 basis=C0_AffineLinearOnCubeWithNodalBasis
         elementQuadrature = CubeGaussQuadrature(nd,3)
         elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,3)
    else:
    	 basis=C0_AffineLinearOnSimplexWithNodalBasis
         elementQuadrature = SimplexGaussQuadrature(nd,3)
         elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
         #elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
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


ns_forceStrongDirichlet = False
backgroundDiffusionFactor = 0.01
sc = 0.5
sc_beta = 1.5
if useMetrics:
    ns_shockCapturingFactor  = sc
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = sc
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = sc_beta
    vof_shockCapturingFactor = sc
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = sc_beta
    rd_shockCapturingFactor  =sc
    rd_lag_shockCapturing = False
    epsFact_density    = 3.
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.
    redist_Newton = True#False
    kappa_shockCapturingFactor = sc
    kappa_lag_shockCapturing = False#True
    kappa_sc_uref = 1.0
    kappa_sc_beta = sc_beta
    dissipation_shockCapturingFactor = sc
    dissipation_lag_shockCapturing = False#True
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = sc_beta
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
    redist_Newton = False#True
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

tolfac = 0.001
ns_nl_atol_res = max(1.0e-6,tolfac*he**2)
vof_nl_atol_res = max(1.0e-6,tolfac*he**2)
ls_nl_atol_res = max(1.0e-6,tolfac*he**2)
mcorr_nl_atol_res = max(1.0e-6,0.1*tolfac*he**2)
rd_nl_atol_res = max(1.0e-6,tolfac*he)
kappa_nl_atol_res = max(1.0e-6,tolfac*he**2)
dissipation_nl_atol_res = max(1.0e-6,tolfac*he**2)
mesh_nl_atol_res = 1.0e-10#max(1.0e-6,tolfac*he**2)
am_nl_atol_res = max(1.0e-6,tolfac*he**2)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure == 4


def twpflowPressure_init(x, t):
    p_L = 0.0
    phi_L = tank_dim[nd-1] - water_level
    phi = x[nd-1] - water_level
    return p_L -g[nd-1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*he,phi)))
