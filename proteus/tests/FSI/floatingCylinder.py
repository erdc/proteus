import numpy as np
from proteus import Domain, Context, Comm
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus import WaveTools as wt
from proteus.Profiling import logEvent
from proteus.mbd import CouplingFSI as fsi
import os
import pychrono

rho_0 = 998.2
nu_0 = 1.004e-6
rho_1 = 1.205
nu_1 = 1.5e-5
sigma_01 = 0.
he = 0.05
tank_dim = [1., 1.]
water_level = 0.5
genMesh = False
rhor = 0.5


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
tank = st.Tank2D(domain, tank_dim)

# CAISSON
radius = 0.1
#caisson = st.Circle(domain,
#                    radius=radius,
#                    coords=(tank_dim[0]/2., water_level+radius-2.*radius*rhor+radius/100.),
#                    barycenter=(tank_dim[0]/2., water_level+radius-2.*radius*rhor+radius/100.),
#                    nPoints=int(np.pi*radius/he))
caisson = st.Rectangle(domain,
                       dim=[2*radius, 2*radius],
                       coords=(tank_dim[0]/2., water_level+radius/10.),
                       barycenter=(tank_dim[0]/2., water_level+radius/10.))
caisson.setHoles([caisson.barycenter[:2]])
caisson.holes_ind = np.array([0])
# let gmsh know that the caisson is IN the tank
tank.setChildShape(caisson, 0)


#  ____                        _                   ____                _ _ _   _
# | __ )  ___  _   _ _ __   __| | __ _ _ __ _   _ / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
# |  _ \ / _ \| | | | '_ \ / _` |/ _` | '__| | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
# | |_) | (_) | |_| | | | | (_| | (_| | |  | |_| | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |____/ \___/ \__,_|_| |_|\__,_|\__,_|_|   \__, |\____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
#                                           |___/
# Boundary Conditions

tank.BC['y+'].setAtmosphere()
tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setFreeSlip()
tank.BC['x-'].setFreeSlip()
tank.BC['sponge'].setNonMaterial()

tank.BC['x-'].setFixedNodes()
tank.BC['x+'].setFixedNodes()
tank.BC['sponge'].setFixedNodes()
tank.BC['y+'].setFixedNodes()  # sliding mesh nodes
tank.BC['y-'].setFixedNodes()  #sliding mesh nodes

for bc in caisson.BC_list:
    bc.setNoSlip()


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
        p_L = 0.0
        phi_L = tank_dim[nd-1] - water_level
        phi = x[nd-1] - water_level
        p = p_L -g[nd-1]*(rho_0*(phi_L - phi)
                          +(rho_1 -rho_0)*(smoothedHeaviside_integral(smoothing,phi_L)
                                                -smoothedHeaviside_integral(smoothing,phi)))
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
class VF_IC:
    def uOfXT(self, x, t):
        return smoothedHeaviside(smoothing,x[nd-1]-water_level)
class PHI_IC:
    def uOfXT(self, x, t):
        return x[nd-1] - water_level

# instanciating the classes for *_p.py files
initialConditions = {'pressure': P_IC(),
                     'vel_u': U_IC(),
                     'vel_v': V_IC(),
                     'vel_w': W_IC(),
                     'vof': VF_IC(),
                     'ncls': PHI_IC(),
                     'rdls': PHI_IC()}

#   ____ _
#  / ___| |__  _ __ ___  _ __   ___
# | |   | '_ \| '__/ _ \| '_ \ / _ \
# | |___| | | | | | (_) | | | | (_) |
#  \____|_| |_|_|  \___/|_| |_|\___/
# Chrono

# System
g = np.array([0., -9.81, 0.])
system = fsi.ProtChSystem()
system.ChSystem.Set_G_acc(pychrono.ChVectorD(g[0], g[1], g[2]))
system.setTimeStep(1e-5)
#system.setCouplingScheme("CSS", prediction="backwardEuler")
# Body
body = fsi.ProtChBody(system=system)
body.attachShape(caisson)
#body.Aij_factor = 1/width
chbod = body.ChBody
x, y, z = caisson.barycenter
pos = pychrono.ChVectorD(x, y, z)
mass = (2.*radius)**2*rho_0*rhor
inertia = pychrono.ChVectorD(1., 1., 1.)
chbod.SetPos(pos)
chbod.SetMass(mass)
chbod.SetInertiaXX(inertia)
#chbod.SetBodyFixed(True)
body.setConstraints(free_x=np.array([0.,1.,0.]), free_r=np.array([0.,0.,0.]))

# body.setInitialRot(rotation_init)
# body.rotation_init=np.array([np.cos(ang/2.), 0., 0., np.sin(ang/2.)*1.])
body.setRecordValues(all_values=True)

#  __  __           _        ___        _   _
# |  \/  | ___  ___| |__    / _ \ _ __ | |_(_) ___  _ __  ___
# | |\/| |/ _ \/ __| '_ \  | | | | '_ \| __| |/ _ \| '_ \/ __|
# | |  | |  __/\__ \ | | | | |_| | |_) | |_| | (_) | | | \__ \
# |_|  |_|\___||___/_| |_|  \___/| .__/ \__|_|\___/|_| |_|___/
#                                |_|


domain.MeshOptions.use_gmsh = genMesh
domain.MeshOptions.genMesh = genMesh
he = he
domain.MeshOptions.he = he
modulepath = os.path.dirname(os.path.abspath(__file__))
mesh_fileprefix=modulepath+'/meshFloatingCylinder'
domain.MeshOptions.setOutputFiles(mesh_fileprefix)
st.assembleDomain(domain)
domain.use_gmsh = False
domain.geofile = mesh_fileprefix



#  _   _                           _
# | \ | |_   _ _ __ ___   ___ _ __(_) ___ ___
# |  \| | | | | '_ ` _ \ / _ \ '__| |/ __/ __|
# | |\  | |_| | | | | | |  __/ |  | | (__\__ \
# |_| \_|\__,_|_| |_| |_|\___|_|  |_|\___|___/
# Numerics

outputStepping = TpFlow.OutputStepping(
    final_time=0.1,
    dt_init=0.01,
    dt_output=0.1,
    nDTout=None,
    dt_fixed=0.01,
)

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(
    ns_model=None,
    ls_model=None,
    nd=domain.nd,
    cfl=0.9,
    outputStepping=outputStepping,
    structured=False,
    he=he,
    nnx=None,
    nny=None,
    nnz=None,
    domain=domain,
    initialConditions=initialConditions,
    boundaryConditions=None, # set with SpatialTools,
    useSuperlu=False,
)

# line below needed for relaxation zones
# (!) hack
m = myTpFlowProblem.Parameters.Models
m.rans2p.auxiliaryVariables += domain.auxiliaryVariables['twp']
myTpFlowProblem.archiveAllSteps = True

myTpFlowProblem.movingDomain = True

params = myTpFlowProblem.Parameters

# MESH PARAMETERS
params.mesh.genMesh = genMesh
params.mesh.he = he

# PHYSICAL PARAMETERS
params.physical.densityA = rho_0  # water
params.physical.densityB = rho_1  # air
params.physical.kinematicViscosityA = nu_0  # water
params.physical.kinematicViscosityB = nu_1  # air
params.physical.gravity = np.array(g)
params.physical.surf_tension_coeff = sigma_01

# MODEL PARAMETERS
ind = -1
m.moveMeshElastic.index = ind+1
ind += 1
m.rans2p.index = ind+1
ind += 1
m.vof.index = ind+1
ind += 1
m.ncls.index = ind+1
ind += 1
m.rdls.index = ind+1
ind += 1
m.mcorr.index = ind+1
ind += 1
m.addedMass.index = ind+1
ind += 1

m.rans2p.auxiliaryVariables += [system]
m.rans2p.p.coefficients.eb_bc_penalty_constant = 10.#/nu_0#Re
m.addedMass.auxiliaryVariables += [system.ProtChAddedMass]
max_flag = 0
max_flag = max(domain.vertexFlags)
max_flag = max(domain.segmentFlags+[max_flag])
max_flag = max(domain.facetFlags+[max_flag])
flags_rigidbody = np.zeros(max_flag+1, dtype='int32')
for s in system.subcomponents:
    if type(s) is fsi.ProtChBody:
        for i in s.boundaryFlags:
            flags_rigidbody[i] = 1
m.addedMass.p.coefficients.flags_rigidbody = flags_rigidbody

