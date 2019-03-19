"""
Non linear waves
"""
from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import WaveTools as wt
import math
import numpy as np
from proteus.ctransportCoefficients import smoothedHeaviside_integral
# For TwoPhaseFlow #
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow

opts=Context.Options([
    # PARAMETERS FOR TwoPhaseFlow #
    ("ns_model",0,"0: rans2p, 1: rans3p"),
    ('ls_model',1,"ls_model = {ncls,clsvof}"),
    ("final_time",8.0,"Final time"),
    ("dt_output",0.1,"How ofter to output solution"),
    ("cfl", 0.5 , "Target cfl"),    
    ("useSuperlu",True,"Use (or not) super lu"),
    # END OF PARAMETERS FOR TwoPhaseFlow #
    # predefined test cases    
    ("water_level", 1.0, "Water level from y=0"),
    # tank
    ("tank_dim", (10.0, 2.0), "Dimensions of the operational domain of tank (l x h)"),
    ("generation", True, "Generate waves at the left boundary (True/False)"),
    ("absorption", True, "Absorb waves at the right boundary (True/False)"),
    ("tank_sponge", (2.0,2.0), "Length of relaxation zones zones in m (left, right)"),
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
    ("gauge_dx", 0.25, "Horizontal spacing of point gauges/column gauges in m"),
    # mesh refinement
    ("refinement", False, "Gradual refinement"),
    ("he", 0.02, "Set characteristic element size in m"),
    ("he_max", 10, "Set maximum characteristic element size in m"),
    ("he_max_water", 10, "Set maximum characteristic in water phase in m"),
    ("refinement_freesurface", 0.05,"Set area of constant refinement around free surface (+/- value) in m"),
    ("refinement_grading", np.sqrt(1.1*4./np.sqrt(3.))/np.sqrt(1.*4./np.sqrt(3)), "Grading of refinement/coarsening (default: 10% volume)"),
    ("gen_mesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("use_gmsh", False, "True: use Gmsh. False: use Triangle/Tetgen")
    ])

# ----- CONTEXT ------ #
# waves
omega = 1.
if opts.waves is True:
    height = opts.wave_height
    mwl = depth = opts.water_level
    direction = opts.wave_dir
    wave = wt.SolitaryWave(waveHeight = height, 
			   mwl = mwl, 
			   depth = depth,
               		   g = np.array(opts.g), 
			   waveDir = direction,
			   trans = np.array([-8*opts.tank_sponge[0], 0., 0.]),
                       	   fast = opts.fast)

# tank options
waterLevel = opts.water_level
tank_dim = opts.tank_dim
tank_sponge = opts.tank_sponge

# ----- DOMAIN ----- #
nd = 2
domain = Domain.PlanarStraightLineGraphDomain()

# refinement
smoothing = opts.he*3

# ----- TANK ------ #
tank = st.Tank2D(domain, tank_dim)

# ----- GENERATION / ABSORPTION LAYERS ----- #
tank.setSponge(x_n=tank_sponge[0], x_p=tank_sponge[1])
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
tank.BC['sponge'].setNonMaterial()
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
    tank.MeshOptions.refineBox(he2, he_max, -tank_sponge[0], tank_dim[0]+tank_sponge[1], waterLevel-box, waterLevel+box)
    tank.MeshOptions.setRefinementFunction(mesh_grading(start='y-{0}'.format(waterLevel-box), he=he2, grading=grading))
    tank.MeshOptions.setRefinementFunction(mesh_grading(start='y-{0}'.format(waterLevel+box), he=he2, grading=grading))
    domain.MeshOptions.LcMax = he_max #coarse grid
    if opts.use_gmsh is True and opts.refinement is True:
        domain.MeshOptions.he = he_max #coarse grid
    else:
        domain.MeshOptions.he = he2 #coarse grid
        domain.MeshOptions.LcMax = he2 #coarse grid
    tank.MeshOptions.refineBox(opts.he_max_water, he_max, -tank_sponge[0], tank_dim[0]+tank_sponge[1], 0., waterLevel)
else:
    domain.MeshOptions.LcMax = opts.he
mr._assembleRefinementOptions(domain)
from proteus import Comm
comm = Comm.get()
if domain.use_gmsh is True:
    mr.writeGeo(domain, 'mesh', append=False)
# ----- END OF REFINEMENT OPTIONS ----- #

#################### 
# FOR TwoPhaseFlow # #Added by mql
####################
# INITIAL CONDITIONS #
class pressure_init_cond:
    def uOfXT(self, x, t):
        p_L = 0.0
        phi_L = tank_dim[nd-1] - waterLevel
        phi = x[nd-1] - waterLevel
        physical_parameters = TpFlow.default_physical_parameters
        clsvof_parameters = TpFlow.default_clsvof_parameters
        g = physical_parameters['gravity']
        rho_0 = physical_parameters['densityA']
        rho_1 = physical_parameters['densityB']
        eps = clsvof_parameters['epsFactHeaviside']*opts.he        
        return p_L -g[nd-1]*(rho_0*(phi_L - phi)+
                             (rho_1 -rho_0)*(smoothedHeaviside_integral(eps,phi_L)
                                             -smoothedHeaviside_integral(eps,phi)))
class zero:
    def uOfXT(self, x, t):
        return 0.0
class clsvof_init_cond:
    def uOfXT(self, x, t):
        return x[nd-1] - waterLevel

############################################
# ***** Create myTwoPhaseFlowProblem ***** #
############################################
outputStepping = TpFlow.OutputStepping(opts.final_time,dt_output=opts.dt_output)
initialConditions = {'pressure': pressure_init_cond(),
                     'pressure_increment': zero(),
                     'vel_u': zero(),
                     'vel_v': zero(),
                     'clsvof': clsvof_init_cond()}
boundaryConditions = {
    # DIRICHLET BCs #
    'pressure_DBC': lambda x, flag: domain.bc[flag].p_dirichlet.init_cython(),
    'pressure_increment_DBC': lambda x, flag: domain.bc[flag].pInc_dirichlet.init_cython(),
    'vel_u_DBC': lambda x, flag: domain.bc[flag].u_dirichlet.init_cython(),
    'vel_v_DBC': lambda x, flag: domain.bc[flag].v_dirichlet.init_cython(),
    'clsvof_DBC': lambda x, flag: domain.bc[flag].vof_dirichlet.init_cython(),
    # ADVECTIVE FLUX BCs #
    'pressure_AFBC': lambda x, flag: domain.bc[flag].p_advective.init_cython(),
    'pressure_increment_AFBC': lambda x, flag: domain.bc[flag].pInc_advective.init_cython(),
    'vel_u_AFBC': lambda x, flag: domain.bc[flag].u_advective.init_cython(),
    'vel_v_AFBC': lambda x, flag: domain.bc[flag].v_advective.init_cython(),
    'clsvof_AFBC': lambda x, flag: domain.bc[flag].vof_advective.init_cython(),
    # DIFFUSIVE FLUX BCs #
    'pressure_increment_DFBC': lambda x, flag: domain.bc[flag].pInc_diffusive.init_cython(),
    'vel_u_DFBC': lambda x, flag: domain.bc[flag].u_diffusive.init_cython(),
    'vel_v_DFBC': lambda x, flag: domain.bc[flag].v_diffusive.init_cython(),
    'clsvof_DFBC': lambda x, flag: None}
myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(ns_model=opts.ns_model,
                                             ls_model=opts.ls_model,
                                             nd=2,
                                             cfl=opts.cfl,
                                             outputStepping=outputStepping,
                                             structured=False,
                                             he=opts.he,
                                             nnx=None,
                                             nny=None,
                                             nnz=None,
                                             domain=domain,
                                             initialConditions=initialConditions,
                                             boundaryConditions=boundaryConditions)

params = myTpFlowProblem.Parameters

# MESH PARAMETERS
params.mesh.genMesh = opts.gen_mesh
params.mesh.he = opts.he

# if opts.ns_model == 0:
#     myTpFlowProblem.Parameters.Models.rans2p['index'] = 0
#     myTpFlowProblem.Parameters.Models.clsvof['index'] = 1
# elif opts.ns_model == 1:
#     myTpFlowProblem.Parameters.Models.clsvof['index'] = 0
#     myTpFlowProblem.Parameters.Models.rans3p['index'] = 1
#     myTpFlowProblem.Parameters.Models.pressureIncrement['index'] = 2
#     myTpFlowProblem.Parameters.Models.pressure['index'] = 3
#     myTpFlowProblem.Parameters.Models.pressureInitial['index'] = 4
