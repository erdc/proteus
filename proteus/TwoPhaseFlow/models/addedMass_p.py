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

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
flags_rigidbody = params.addedMass['flags_rigidbody']

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
V_model = params.rans2p['index'] or params.rans3p['index']

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
dirichletConditions = {0: lambda x, flag: None}
advectiveFluxBoundaryConditions = {}
def getFlux_am(x, flag):
    #the unit rigid motions will applied internally
    #leave this set to zero
    return lambda x,t: 0.0
diffusiveFluxBoundaryConditions = {0: {0: getFlux_am}}

