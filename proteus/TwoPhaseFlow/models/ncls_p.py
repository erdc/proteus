from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus.mprans import NCLS
from proteus import Context

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()

myTpFlowProblem = ct.myTpFlowProblem 
params = myTpFlowProblem.Parameters
initialConditions   = myTpFlowProblem.initialConditions
boundaryConditions  = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd

# DOMAIN #
domain = myTpFlowProblem.domain

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
useMetrics = params.ncls['useMetrics']
checkMass = params.ncls['checkMass']
sc_uref = params.ncls['sc_uref']
sc_beta = params.ncls['sc_beta']
epsFact = params.ncls['epsFact']

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
ME_model = params.ncls['index']
assert ME_model != None, 'ls model index was not set!'
if params.rans2p['index'] is not None:
    V_model = params.rans2p['index']
elif params.rans3p['index'] is not None:
    V_model = params.rans3p['index']
else:
    assert params.rans2p['index'] is not None or params.rans3p['index'] is not None, 'RANS2P or RANS3P must be used with VOF'
RD_model = params.rdls['index']

LevelModelType = NCLS.LevelModel
coefficients = NCLS.Coefficients(V_model=V_model,
                                 RD_model=RD_model,
                                 ME_model=ME_model,
                                 checkMass=checkMass,
                                 useMetrics=useMetrics,
                                 epsFact=epsFact,
                                 sc_uref=sc_uref,
                                 sc_beta=sc_beta,
                                 movingDomain=ct.movingDomain)

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions = {0: initialConditions['ncls']}

# ***************************************** #
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
dirichletConditions = {0: lambda x, flag: None}
advectiveFluxBoundaryConditions = {}
diffusiveFluxBoundaryConditions = {0: {}}
