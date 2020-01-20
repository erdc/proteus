import numpy as np
from proteus import Domain
from proteus import Context
from proteus.mprans import SpatialTools as st
from proteus.mbd import CouplingFSI as fsi
import pychrono
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow

addedMass = False

#   ____            _            _      ___        _   _
#  / ___|___  _ __ | |_ _____  _| |_   / _ \ _ __ | |_(_) ___  _ __  ___
# | |   / _ \| '_ \| __/ _ \ \/ / __| | | | | '_ \| __| |/ _ \| '_ \/ __|
# | |__| (_) | | | | ||  __/>  <| |_  | |_| | |_) | |_| | (_) | | | \__ \
#  \____\___/|_| |_|\__\___/_/\_\\__|  \___/| .__/ \__|_|\___/|_| |_|___/
#                                           |_|
# Context options
# used in command line directly with option -C
# e.g.: parun [...] -C "g=(0.,-9.81,0.) rho_0=998.2 genMesh=False"
#
# only change/add context options in the "other options section"
# other sections have variables used in _p and _n files

context_options = []
# physical constants
context_options += [
    ("densityA", 998.2, "Water density"),
    ("densityB", 1.004e-6, "Water kinematic viscosity m/sec^2"),
    ("kinematicViscosityA", 1.205, "Air Densiy"),
    ("kinematicViscosityB", 1.5e-5, "Air kinematic viscosity m/sec^2"),
    ("surf_tension_coeff", 0., "Surface tension"),
    ("gravity", (0, -9.81, 0.), "Gravitational acceleration vector"),
    ]
# run time options
context_options += [
    ("T", 1., "Simulation time in s"),
    ("dt_init", 0.001 ,"Value of initial time step"),
    ("dt_fixed", 0.5, "Value of maximum time step"),
    ("archiveAllSteps", False, "archive every steps"),
    ("dt_output", 1., "number of saves per second"),
    ("cfl", 0.9 , "Target CFL value"),
    ]
# mesh options
context_options += [
    ("he", 1. , "Characteristic element size"),
    ]
# other options
context_options += [
    ("radius", 0.5, "Radius of cylinder"),
    ("use_ball_as_particle", 0., "use ball method"),
    ]

opts = Context.Options(context_options)
he = opts.he
radius = opts.radius

#  ____                        _
# |  _ \  ___  _ __ ___   __ _(_)_ __
# | | | |/ _ \| '_ ` _ \ / _` | | '_ \
# | |_| | (_) | | | | | | (_| | | | | |
# |____/ \___/|_| |_| |_|\__,_|_|_| |_|
# Domain
# All geometrical options go here (but not mesh options)

domain = Domain.PlanarStraightLineGraphDomain()

# ----- SHAPES ----- #

# TANK
tank_dim = np.array([3., 10.])
tank_x = tank_dim[0]
tank_y = tank_dim[1]
water_level=tank_y*2
tank = st.Tank2D(domain, [tank_x, tank_y])

#  ____                        _                   ____                _ _ _   _
# | __ )  ___  _   _ _ __   __| | __ _ _ __ _   _ / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
# |  _ \ / _ \| | | | '_ \ / _` |/ _` | '__| | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
# | |_) | (_) | |_| | | | | (_| | (_| | |  | |_| | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |____/ \___/ \__,_|_| |_|\__,_|\__,_|_|   \__, |\____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
#                                           |___/
# Boundary Conditions


#tank.BC['y+'].setFreeSlip()
tank.BC['y+'].setAtmosphere()
#tank.BC['y-'].setFreeSlip()
#tank.BC['x+'].setFreeSlip()
#tank.BC['x-'].setFreeSlip()
#tank.BC['y+'].setNoSlip()
tank.BC['y-'].setNoSlip()
tank.BC['x+'].setNoSlip()
tank.BC['x-'].setNoSlip()
#tank.BC['sponge'].setNonMaterial()


#   ____ _
#  / ___| |__  _ __ ___  _ __   ___
# | |   | '_ \| '__/ _ \| '_ \ / _ \
# | |___| | | | | | (_) | | | | (_) |
#  \____|_| |_|_|  \___/|_| |_|\___/
# Chrono
g = np.array(opts.gravity)

# System
system = fsi.ProtChSystem()
system.ChSystem.Set_G_acc(pychrono.ChVectorD(g[0], g[1], g[2]))
system.setTimeStep(1e-5)
# Body
body = fsi.ProtChBody(system=system)
body.setBoundaryFlags([0])  # index of particle
#
chbod = body.ChBody
pos = pychrono.ChVectorD(0.5*tank_x, 0.7*tank_y, 0.)
rot = pychrono.ChQuaternionD(1., 0., 0., 0.)
mass = 2000.0#3.14*rho_0*1.5
inertia = pychrono.ChVectorD(1., 1., 1.)
chbod.SetPos(pos)
chbod.SetRot(rot)
chbod.SetMass(mass)
chbod.SetInertiaXX(inertia)
# chbod.SetBodyFixed(True)
body.setConstraints(free_x=np.array([0., 1., 0.]), free_r=np.array([0., 0., 0.]))
body.setRecordValues(all_values=True)


def sdf(t, x):
    dist = np.sqrt((x[0]**2+x[1]**2+x[2]**2))
    if dist < radius:
        return -dist, (0., -1., 0.)
    else:
        return dist, (0., 1., 0.)

body.setIBM(True,
            radiusIBM=radius,  # used only when particle are balls
            sdfIBM=sdf,  # used only when particle not used as balls
)


#  ___       _ _   _       _    ____                _ _ _   _
# |_ _|_ __ (_) |_(_) __ _| |  / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
#  | || '_ \| | __| |/ _` | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
#  | || | | | | |_| | (_| | | | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |___|_| |_|_|\__|_|\__,_|_|  \____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
# Initial Conditions

from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
smoothing = 1.5 * he
nd = domain.nd

class P_IC:
    def uOfXT(self, x, t):
        p = -opts.densityA*g[nd-1]*(tank_dim[1] - x[1])
        return p
class U_IC:
    def uOfXT(self, x, t):
        return 0.0
class V_IC:
    def uOfXT(self, x, t):
        return 0.0
class W_IC:
    def uOfXT(self, x, t):
        return 0.0
# instanciating the classes for *_p.py files
initialConditions = {'pressure': P_IC(),
                     'vel_u': U_IC(),
                     'vel_v': V_IC(),
                     'vel_w': W_IC()}

#  __  __           _        ___        _   _
# |  \/  | ___  ___| |__    / _ \ _ __ | |_(_) ___  _ __  ___
# | |\/| |/ _ \/ __| '_ \  | | | | '_ \| __| |/ _ \| '_ \/ __|
# | |  | |  __/\__ \ | | | | |_| | |_) | |_| | (_) | | | \__ \
# |_|  |_|\___||___/_| |_|  \___/| .__/ \__|_|\___/|_| |_|___/
#                                |_|


domain.MeshOptions.use_gmsh = False
domain.MeshOptions.genMesh = False
domain.MeshOptions.he = he
mesh_fileprefix='mesh'+str(int(1000*he))
domain.MeshOptions.setOutputFiles(mesh_fileprefix)
st.assembleDomain(domain)
domain.use_gmsh = False
domain.geofile = mesh_fileprefix

domain.BCbyFlag[0].u_diffusive.uOfXT = lambda x, t: 0.
domain.BCbyFlag[0].v_diffusive.uOfXT = lambda x, t: 0.

#  _   _                           _
# | \ | |_   _ _ __ ___   ___ _ __(_) ___ ___
# |  \| | | | | '_ ` _ \ / _ \ '__| |/ __/ __|
# | |\  | |_| | | | | | |  __/ |  | | (__\__ \
# |_| \_|\__,_|_| |_| |_|\___|_|  |_|\___|___/
# Numerics

outputStepping = TpFlow.OutputStepping(
    final_time=opts.T,
    dt_init=opts.dt_init,
    dt_output=opts.dt_output,
    nDTout=None,
    dt_fixed=opts.dt_fixed,
)

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(
    ns_model=None,
    ls_model=None,
    nd=domain.nd,
    cfl=opts.cfl,
    outputStepping=outputStepping,
    structured=False,
    he=he,
    nnx=None,
    nny=None,
    nnz=None,
    domain=domain,
    initialConditions=initialConditions,
    boundaryConditions=None, # set with SpatialTools,
    useSuperlu=True,
)

# line below needed for relaxation zones
# (!) hack
m = myTpFlowProblem.Parameters.Models
m.rans2p.auxiliaryVariables += domain.auxiliaryVariables['twp']
myTpFlowProblem.archiveAllSteps = opts.archiveAllSteps

myTpFlowProblem.movingDomain = False

params = myTpFlowProblem.Parameters

# MESH PARAMETERS
params.mesh.genMesh = True
params.mesh.he = he

# PHYSICAL PARAMETERS

params.physical.densityA = opts.densityA  # water
params.physical.densityB = opts.densityB  # air
params.physical.kinematicViscosityA = opts.kinematicViscosityA  # water
params.physical.kinematicViscosityB = opts.kinematicViscosityB  # air
params.physical.gravity = np.array(opts.gravity)
params.physical.surf_tension_coeff = 0.

# MODEL PARAMETERS
ind = -1
m.rans2p.index = ind+1
ind += 1
if addedMass is True:
    m.addedMass.index = ind+1
    ind += 1
m.rans2p.auxiliaryVariables += [system]
# m.rans2p.p.CoefficientsOptions.particle_epsFact = 1.5
m.rans2p.p.CoefficientsOptions.use_ball_as_particle = opts.use_ball_as_particle
m.rans2p.p.CoefficientsOptions.nParticles = 1
#m.rans2p.p.CoefficientsOptions.particle_alpha = 1000.
#m.rans2p.p.CoefficientsOptions.particle_beta = 1000.
#m.rans2p.p.CoefficientsOptions.particle_penalty_constant = 100.
#m.rans2p.p.CoefficientsOptions.particle_sdfList = [lambda t, x: body.getDynamicSDF(t, x)]
#m.rans2p.p.CoefficientsOptions.particle_velocityList = [lambda t, x: body.getVelocity()]
# m.rans2p.p.CoefficientsOptions.useExact = False
# m.rans2p.p.CoefficientsOptions.ghost_penalty_constant = 0.1
#m.rans2p.p.CoefficientsOptions.ball_radius = np.array([radius], 'd')
#m.rans2p.p.CoefficientsOptions.ball_velocity = np.zeros((1, 3), 'd')
#m.rans2p.p.CoefficientsOptions.ball_angular_velocity = np.zeros((1, 3), 'd')
#m.rans2p.p.CoefficientsOptions.ball_center = np.array([[pos.x, pos.y, pos.z]], 'd')

# m.rans2p.p.CoefficientsOptions.MOMENTUM_SGE = 1.
# m.rans2p.p.CoefficientsOptions.PRESSURE_SGE = 1.
# m.rans2p.p.CoefficientsOptions.VELOCITY_SGE = 1.
