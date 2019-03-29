from __future__ import print_function
from __future__ import division
from past.utils import old_div
import numpy as np
from proteus import Domain
from proteus.mprans import SpatialTools as st
from proteus.mbd import CouplingFSI as fsi
from proteus.TwoPhaseFlow import TwoPhaseFlowProblem as tpf
import pychrono as chrono

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

system = fsi.ProtChSystem()
system.ChSystem.Set_G_acc(chrono.ChVectorD(g[0], g[1], 0.))
body = fsi.ProtChBody(system=system)
body.attachShape(rect)
body.ChBody.SetMass(500.)
body.ChBody.SetBodyFixed(True)  # fixing body

# OTHER PARAMS
st.assembleDomain(domain)

#  ___       _ _   _       _    ____                _ _ _   _
# |_ _|_ __ (_) |_(_) __ _| |  / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
#  | || '_ \| | __| |/ _` | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
#  | || | | | | |_| | (_| | | | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |___|_| |_|_|\__|_|\__,_|_|  \____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
# Initial Conditions

nd = domain.nd

class PerturbedSurface_p:
    def uOfXT(self, x, t):
        p_L = 0.0
        phi = x[nd-1] - tank_dim[nd-1]
        return p_L-g[nd-1]*(rho_0*phi)

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0


# instanciating the classes for *_p.py files
initialConditions = {'pressure': PerturbedSurface_p(),
                     'vel_u': AtRest(),
                     'vel_v': AtRest()}


#  _   _                           _
# | \ | |_   _ _ __ ___   ___ _ __(_) ___ ___
# |  \| | | | | '_ ` _ \ / _ \ '__| |/ __/ __|
# | |\  | |_| | | | | | |  __/ |  | | (__\__ \
# |_| \_|\__,_|_| |_| |_|\___|_|  |_|\___|___/
# Numerics


outputStepping = tpf.OutputStepping(
    final_time=0.002,
    dt_init=0.001,
    dt_output=0.001,
    nDTout=None,
    dt_fixed=0.001,
)

myTpFlowProblem = tpf.TwoPhaseFlowProblem(
    ns_model=None,
    ls_model=None,
    nd=domain.nd,
    cfl=0.4,
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
myTpFlowProblem.movingDomain = False

params = myTpFlowProblem.Parameters

# PHYSICAL PARAMETERS
params.physical.densityA = rho_0  # water
params.physical.densityB = rho_1  # air
params.physical.kinematicViscosityA = nu_0  # water
params.physical.kinematicViscosityB = nu_1  # air
params.physical.surf_tension_coeff = sigma_01

# MODEL PARAMETERS
m = params.Models
m.rans2p.index = 0
m.addedMass.index = 1

# auxiliary variables
m.rans2p.auxiliaryVariables += [system]
m.addedMass.auxiliaryVariables += [system.ProtChAddedMass]

flags_rigidbody = np.zeros(20)
for key in rect.boundaryTags_global:
    flags_rigidbody[rect.boundaryTags_global[key]] = 1.

max_flag = 0
max_flag = max(domain.vertexFlags)
max_flag = max(domain.segmentFlags+[max_flag])
max_flag = max(domain.facetFlags+[max_flag])
flags_rigidbody = np.zeros(max_flag+1, dtype='int32')
for s in system.subcomponents:
    if type(s) is fsi.ProtChBody:
        for i in range(s.i_start, s.i_end):
            flags_rigidbody[i] = 1
m.addedMass.p.CoefficientsOptions.flags_rigidbody = flags_rigidbody
