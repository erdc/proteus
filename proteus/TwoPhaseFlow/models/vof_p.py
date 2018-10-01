from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus.mprans import VOF
from proteus import Context

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()

myTpFlowProblem = ct.myTpFlowProblem 
initialConditions = myTpFlowProblem.initialConditions
boundaryConditions = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd
movingDomain = myTpFlowProblem.movingDomain

# DOMAIN #
domain = myTpFlowProblem.domain

params = myTpFlowProblem.Parameters
mparams = params.Models # model parameters
myparams = params.Models.vof
pparams = params.physical # physical parameters

# MESH #
meshparams = params.mesh
genMesh = meshparams.genMesh

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
ME_model = mparams.vof.index
assert ME_model is not None, 'vof model index was not set!'
LS_model = mparams.ncls.index
if mparams.rans2p.index is not None:
    V_model = mparams.rans2p.index
elif mparams.rans3p.index is not None:
    V_model = mparams.rans3p.index
else:
    assert mparams.rans2p.index is not None or params.rans3p.index is not None, 'RANS2P or RANS3P must be used with VOF'
RD_model = mparams.rdls.index

LevelModelType = VOF.LevelModel
coefficients = VOF.Coefficients(LS_model=LS_model,
                                V_model=V_model,
                                RD_model=RD_model,
                                ME_model=ME_model,
                                checkMass=myparams.checkMass,
                                useMetrics=myparams.useMetrics,
                                epsFact=myparams.epsFact,
                                sc_uref=myparams.sc_uref,
                                sc_beta=myparams.sc_beta,
                                movingDomain=movingDomain)

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions = {0: initialConditions['vof']}

# ***************************************** #
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
if domain.useSpatialTools is False or myTpFlowProblem.useBoundaryConditionsModule is False:
    dirichletConditions = {0: boundaryConditions['vof_DBC']}
    advectiveFluxBoundaryConditions = {0: boundaryConditions['vof_AFBC']}
else:
    dirichletConditions = {0: lambda x, flag: domain.bc[flag].vof_dirichlet.init_cython()}
    advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].vof_advective.init_cython()}
diffusiveFluxBoundaryConditions = {0: {}}
