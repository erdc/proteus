"""
Non linear waves
"""
from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import WaveTools as wt
import math
import numpy as np

opts=Context.Options([
    # predefined test cases
    ("water_level", 1.0, "Water level from y=0"),
    # tank
    ("tank_dim", (15.0, 2.0), "Dimensions of the operational domain of tank (l x h)"),
    ("generation", False, "Generate waves at the left boundary (True/False)"),
    ("absorption", False, "Absorb waves at the right boundary (True/False)"),
    ("tank_sponge", (0., 0.), "Length of relaxation zones zones in m (left, right)"),
    ("free_slip", True, "Should tank walls have free slip conditions "
                        "(otherwise, no slip conditions will be applied)."),
    # gravity
    ("g", [0, -9.81, 0], "Gravity vector in m/s^2"),
    # waves
    ("waves", True, "Generate waves (True/False)"),
    ("wave_height", 0.45, "Height of the waves in s"),
    ("wave_dir", (1.,0.,0.), "Direction of the waves (from left boundary)"),
    ("wave_type", 'solitaryWave', "type of wave"),
    ("fast", False, "switch for fast cosh calculations in WaveTools"),
    # gauges
    #("gauge_output", True, "Places Gauges in tank (5 per wavelength)"),
    ("point_gauge_output", True, "Produce point gauge output"),
    ("column_gauge_output", True, "Produce column gauge output"),
    ("gauge_dx", 5., "Horizontal spacing of point gauges/column gauges in m"),
    # mesh refinement
    ("refinement", False, "Gradual refinement"),
    ("he", 0.02, "Set characteristic element size in m"),
    ("he_max", 10, "Set maximum characteristic element size in m"),
    ("he_max_water", 10, "Set maximum characteristic in water phase in m"),
    ("refinement_freesurface", 0.05,"Set area of constant refinement around free surface (+/- value) in m"),
    ("refinement_caisson", 0.,"Set area of constant refinement (Box) around potential structure (+/- value) in m"),
    ("refinement_grading", np.sqrt(1.1*4./np.sqrt(3.))/np.sqrt(1.*4./np.sqrt(3)), "Grading of refinement/coarsening (default: 10% volume)"),
    # numerical options
    ("gen_mesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("use_gmsh", False, "True: use Gmsh. False: use Triangle/Tetgen"),
    ("movingDomain", False, "True/False"),
    ("T", 3.0, "Simulation time in s"),
    ("dt_init", 0.001, "Initial time step in s"),
    ("dt_fixed", None, "Fixed (maximum) time step"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ("cfl", 0.5 , "Target cfl"),
    ("nsave",  5, "Number of time steps to save per second"),
    ("useRANS", 0, "RANS model"),
    ])

# ----- CONTEXT ------ #

# waves
omega = 1.
if opts.waves is True:
    height = opts.wave_height
    mwl = depth = opts.water_level
    direction = opts.wave_dir
    wave = wt.SolitaryWave(	waveHeight = height, 
				mwl = mwl, 
				depth = depth,
               			g = np.array(opts.g), 
				waveDir = direction,
				trans = np.array([opts.tank_dim[0]/2, 0., 0.]),
                       		fast = opts.fast
			  )

# tank options
waterLevel = opts.water_level
tank_dim = opts.tank_dim
tank_sponge = opts.tank_sponge

# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()

# refinement
smoothing = opts.he*3

# ----- TANK ------ #

tank = st.Tank2D(domain, tank_dim)

# ----- GENERATION / ABSORPTION LAYERS ----- #

#tank.setSponge(x_n=tank_sponge[0], x_p=tank_sponge[1])
dragAlpha = 5.*omega/1e-6
 
if opts.generation:
    tank.setGenerationZones(x_n=True, waves=wave, dragAlpha=dragAlpha, smoothing = smoothing)
if opts.absorption:
    tank.setAbsorptionZones(x_p=True, dragAlpha = dragAlpha)

# ----- BOUNDARY CONDITIONS ----- #

# Waves
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=smoothing, vert_axis=1)

# open top
tank.BC['y+'].setAtmosphere()

if opts.free_slip:
    tank.BC['y-'].setFreeSlip()
    tank.BC['x+'].setFreeSlip()
    if not opts.generation:
        tank.BC['x-'].setFreeSlip()
else:  # no slip
    tank.BC['y-'].setNoSlip()
    tank.BC['x+'].setNoSlip()

# sponge
#tank.BC['sponge'].setNonMaterial()

for bc in tank.BC_list:
    bc.setFixedNodes()

# ----- GAUGES ----- #

column_gauge_locations = []
point_gauge_locations = []

if opts.point_gauge_output or opts.column_gauge_output:
    gauge_y = waterLevel - 0.5 * depth
    number_of_gauges = tank_dim[0] / opts.gauge_dx + 1
    for gauge_x in np.linspace(0, tank_dim[0], int(number_of_gauges)):
        point_gauge_locations.append((gauge_x, gauge_y, 0), )
        column_gauge_locations.append(((gauge_x, 0., 0.),
                                       (gauge_x, tank_dim[1], 0.)))

if opts.point_gauge_output:
    tank.attachPointGauges('twp',
                           gauges=((('p',), point_gauge_locations),),
                           fileName='pressure_gaugeArray.csv')

if opts.column_gauge_output:
    tank.attachLineIntegralGauges('vof',
                                  gauges=((('vof',), column_gauge_locations),),
                                  fileName='column_gauges.csv')

# ----- ASSEMBLE DOMAIN ----- #

domain.MeshOptions.use_gmsh = opts.use_gmsh
domain.MeshOptions.genMesh = opts.gen_mesh
domain.MeshOptions.he = opts.he
domain.use_gmsh = opts.use_gmsh
st.assembleDomain(domain)

# ----- REFINEMENT OPTIONS ----- #

import MeshRefinement as mr
#domain.MeshOptions = mr.MeshOptions(domain)
tank.MeshOptions = mr.MeshOptions(tank)
if opts.refinement:
    grading = opts.refinement_grading
    he2 = opts.he
    def mesh_grading(start, he, grading):
        return '{0}*{2}^(1+log((-1/{2}*(abs({1})-{0})+abs({1}))/{0})/log({2}))'.format(he, start, grading)
    he_max = opts.he_max
    # he_fs = he2
    ecH = 3.
    if opts.refinement_freesurface > 0:
        box = opts.refinement_freesurface
    else:
        box = ecH*he2
#    tank.MeshOptions.refineBox(he2, he_max, -tank_sponge[0], tank_dim[0]+tank_sponge[1], waterLevel-box, waterLevel+box)
    tank.MeshOptions.setRefinementFunction(mesh_grading(start='y-{0}'.format(waterLevel-box), he=he2, grading=grading))
    tank.MeshOptions.setRefinementFunction(mesh_grading(start='y-{0}'.format(waterLevel+box), he=he2, grading=grading))
    domain.MeshOptions.LcMax = he_max #coarse grid
    if opts.use_gmsh is True and opts.refinement is True:
        domain.MeshOptions.he = he_max #coarse grid
    else:
        domain.MeshOptions.he = he2 #coarse grid
        domain.MeshOptions.LcMax = he2 #coarse grid
#    tank.MeshOptions.refineBox(opts.he_max_water, he_max, -tank_sponge[0], tank_dim[0]+tank_sponge[1], 0., waterLevel)
else:
    domain.MeshOptions.LcMax = opts.he
mr._assembleRefinementOptions(domain)
from proteus import Comm
comm = Comm.get()
if domain.use_gmsh is True:
    mr.writeGeo(domain, 'mesh', append=False)



##########################################
# Numerical Options and other parameters #
##########################################

rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = [0., -9.81]




from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral


#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
movingDomain=opts.movingDomain
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = 10.0/nu_0#Re
dt_init = opts.dt_init
T = opts.T
nDTout = int(opts.T*opts.nsave)
timeIntegration = opts.timeIntegration
if nDTout > 0:
    dt_out= (T-dt_init)/nDTout
else:
    dt_out = 0
runCFL = opts.cfl
dt_fixed = opts.dt_fixed

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
useRANS = opts.useRANS # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega, 1998
            # 3 -- K-Omega, 1988
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


# Numerical parameters
sc = 0.5 # default: 0.5. Test: 0.25
sc_beta = 1.5 # default: 1.5. Test: 1.
epsFact_consrv_diffusion = 1. # default: 1.0. Test: 0.1
ns_forceStrongDirichlet = False
backgroundDiffusionFactor=0.01
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
    epsFact_consrv_diffusion = epsFact_consrv_diffusion
    redist_Newton = True#False
    kappa_shockCapturingFactor = sc
    kappa_lag_shockCapturing = True
    kappa_sc_uref = 1.0
    kappa_sc_beta = sc_beta
    dissipation_shockCapturingFactor = sc
    dissipation_lag_shockCapturing = True
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

ns_nl_atol_res = max(1.0e-6,0.001*domain.MeshOptions.he**2)
vof_nl_atol_res = max(1.0e-6,0.001*domain.MeshOptions.he**2)
ls_nl_atol_res = max(1.0e-6,0.001*domain.MeshOptions.he**2)
mcorr_nl_atol_res = max(1.0e-6,0.0001*domain.MeshOptions.he**2)
rd_nl_atol_res = max(1.0e-6,0.01*domain.MeshOptions.he)
kappa_nl_atol_res = max(1.0e-6,0.001*domain.MeshOptions.he**2)
dissipation_nl_atol_res = max(1.0e-6,0.001*domain.MeshOptions.he**2)
mesh_nl_atol_res = max(1.0e-6,0.001*domain.MeshOptions.he**2)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure == 4

def twpflowPressure_init(x, t):
    p_L = 0.0
    phi_L = tank_dim[nd-1] - waterLevel
    phi = x[nd-1] - waterLevel
    return p_L -g[nd-1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi)))
