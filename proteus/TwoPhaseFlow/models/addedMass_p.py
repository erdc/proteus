from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus import Context
from proteus.mprans import AddedMass

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()

myTpFlowProblem = ct.myTpFlowProblem 
params = myTpFlowProblem.Parameters
initialConditions = myTpFlowProblem.initialConditions
boundaryConditions = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd
movingDomain = myTpFlowProblem.movingDomain

# DOMAIN #
domain = myTpFlowProblem.domain

# MESH #
meshparams = params.mesh
genMesh = meshparams.genMesh

mparams = params.Models
pparams = params.physical
myparams = mparams.addedMass

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
flags_rigidbody = myparams.flags_rigidbody

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
if mparams.rans2p.index is not None:
    V_model = mparams.rans2p.index
elif mparams.rans3p.index is not None:
    V_model = mparams.rans3p.index
else:
    assert mparams.rans2p.index is not None or mparams.rans3p.index is not None, 'RANS2P or RANS3P must be used with addedMass'

LevelModelType = AddedMass.LevelModel

coefficients = AddedMass.Coefficients(nd=nd,
                                      V_model=V_model,
                                      barycenters=domain.barycenters,
                                      flags_rigidbody=flags_rigidbody)

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
class dp_IC:
    def uOfXT(self, x, t):
        return 0.0
initialConditions = {0: dp_IC()}

# ***************************************** #
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
dirichletConditions = {0: lambda x, flag: domain.bc[flag].pAddedMass_dirichlet.init_cython()}
advectiveFluxBoundaryConditions = {}
def getFlux_am(x, flag):
    #the unit rigid motions will applied internally
    #leave this set to zero
    return lambda x,t: 0.0
diffusiveFluxBoundaryConditions = {0: {0: getFlux_am}}

