from proteus.mbd import ChRigidBody as crb
from proteus.mbd import pyChronoCore as pych

import numpy as np

g = np.array([0., -9.81, 0.])

system = crb.ProtChSystem(gravity=g)
system.setTimeStep(0.01)
body = crb.ProtChBody(system=system)

vec = pych.ChVector(0., 2., 0.)
inertia = pych.ChVector(1., 1., 2.)
rot = pych.ChQuaternion(1., 0., 0., 0.)

body.ChBody.SetPos(vec)
body.ChBody.SetRot(rot)
body.ChBody.SetMass(50.)
body.ChBody.SetInertiaXX(inertia)
body.setRecordValues(all_values=True)

body.setConstraints(free_x=np.ones(3), free_r=np.ones(3))

system.calculate_init()
for i in range(10):
    system.calculate(proteus_dt=1.)
