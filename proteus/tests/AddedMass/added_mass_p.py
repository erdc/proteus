from proteus.default_p import *
from proteus.mprans import AddedMass
from proteus import Context
name = "added_mass"
ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions
genMesh = mesh.genMesh
movingDomain = ct.movingDomain
T = ct.T

LevelModelType = AddedMass.LevelModel

coefficients = AddedMass.Coefficients(nd=nd,
                                      V_model=int(ct.movingDomain),
                                      barycenters=domain.barycenters,
                                      flags_rigidbody=ct.flags_rigidbody)

def getDBC_am(x,flag):
    return None

dirichletConditions = {0:getDBC_am}

advectiveFluxBoundaryConditions = {}

def getFlux_am(x, flag):
    #the unit rigid motions will applied internally
    #leave this set to zero
    return lambda x,t: 0.0

diffusiveFluxBoundaryConditions = {0: {0:getFlux_am}}

class dp_IC:
    def uOfXT(self, x, t):
        return 0.0

initialConditions = {0: dp_IC()}
