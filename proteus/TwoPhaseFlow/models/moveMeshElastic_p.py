from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus.mprans import MoveMesh
import numpy as np
from proteus import Context

# *********************************************** #
# ********** READ FROM myTpFlowProblem ********** #
# *********************************************** #
ct = Context.get()

myTpFlowProblem = ct.myTpFlowProblem 
initialConditions   = myTpFlowProblem.initialConditions
boundaryConditions  = myTpFlowProblem.boundaryConditions
nd = myTpFlowProblem.nd
movingDomain = myTpFlowProblem.movingDomain

# DOMAIN #
domain = myTpFlowProblem.domain

params = myTpFlowProblem.Parameters
mparams = params.Models # model parameters
pparams = params.physical # physical parameters

# MESH #
meshparams = params.mesh
genMesh = meshparams.genMesh

# ****************************************** #
# ********** NUMERICAL PARAMETERS ********** #
# ****************************************** #
nMediaTypes = len(domain.regionFlags)  # (!) should be region flags
smTypes     = np.zeros((nMediaTypes+1, 2), 'd')
smFlags     = np.zeros((nMediaTypes+1,), 'i')
smTypes[:, 0] = 1.0    # e
smTypes[:, 1] = 0.3    # nu

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
ME_model = mparams.moveMeshElastic.index
assert ME_model is not None, 'vof model index was not set!'
if mparams.rans2p.index is not None:
    V_model = mparams.rans2p.index
elif mparams.rans3p.index is not None:
    V_model = mparams.rans3p.index
else:
    assert mparams.rans2p.index is not None or mparams.rans3p.index is not None, 'RANS2P or RANS3P must be used with VOF'

LevelModelType = MoveMesh.LevelModel
coefficients = MoveMesh.Coefficients(nd=nd,
                                     V_model=V_model,
                                     modelType_block=smFlags,
                                     modelParams_block=smTypes,
                                     meIndex=ME_model)

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions = None
analyticalSolution = {}

if domain.useSpatialTools is False or myTpFlowProblem.useBoundaryConditionsModule is False:
    dirichletConditions = {0: boundaryConditions['hx'],
                           1: boundaryConditions['hy']}
    stressFluxBoundaryConditions = {0: boundaryConditions['u_stress'],
                                    1: boundaryConditions['v_stress']}
    if nd == 3:
        dirichletConditions[2] = boundaryConditions['hz']
        stressFluxBoundaryConditions[2] = boundaryConditions['w_stress']

else:
    dirichletConditions = {0: lambda x, flag: domain.bc[flag].hx_dirichlet.init_cython(),
                           1: lambda x, flag: domain.bc[flag].hy_dirichlet.init_cython()}
    stressFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].u_stress.init_cython(),
                                    1: lambda x, flag: domain.bc[flag].v_stress.init_cython()}
    if nd == 3:
        dirichletConditions[2] = lambda x, flag: domain.bc[flag].hz_dirichlet.init_cython()
        stressFluxBoundaryConditions[2] = lambda x, flag: domain.bc[flag].w_stress.init_cython()

fluxBoundaryConditions = {0: 'noFlow',
                          1: 'noFlow'}

advectiveFluxBoundaryConditions = {}

diffusiveFluxBoundaryConditions = {0: {},
                                   1: {}}

if nd == 3:
    fluxBoundaryConditions[2] = 'noFlow'
    diffusiveFluxBoundaryConditions[2] = {}

