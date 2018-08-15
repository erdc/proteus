from proteus.default_p import *
from proteus.mprans import MoveMesh
import numpy as np
from proteus import Context

ct = Context.get()
domain = ct.domain
nd = ct.domain.nd
mesh = domain.MeshOptions


genMesh = mesh.genMesh
movingDomain = ct.movingDomain
T = ct.T

initialConditions = None

analyticalSolution = {}

nMediaTypes = len(domain.regionFlags)  # (!) should be region flags
smTypes     = np.zeros((nMediaTypes+1, 2), 'd')
smFlags     = np.zeros((nMediaTypes+1,), 'i')
smTypes[:, 0] = 1.0    # E
smTypes[:, 1] = 0.3    # nu

LevelModelType = MoveMesh.LevelModel
coefficients = MoveMesh.Coefficients(nd=ct.domain.nd,
                                     V_model=1,
                                     modelType_block=smFlags,
                                     modelParams_block=smTypes,
                                     meIndex=0)


dirichletConditions = {0: lambda x, flag: domain.bc[flag].hx_dirichlet.init_cython(),
                       1: lambda x, flag: domain.bc[flag].hy_dirichlet.init_cython()}

fluxBoundaryConditions = {0: 'noFlow',
                          1: 'noFlow'}

advectiveFluxBoundaryConditions = {}

diffusiveFluxBoundaryConditions = {0: {},
                                   1: {}}

stressFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].u_stress,
                                1: lambda x, flag: domain.bc[flag].v_stress}

if nd == 3:
    dirichletConditions[2] = lambda x, flag: domain.bc[flag].hz_dirichlet.init_cython()
    fluxBoundaryConditions[2] = 'noFlow'
    diffusiveFluxBoundaryConditions[2] = {}
    stressFluxBoundaryConditions[2] = lambda x, flag: domain.bc[flag].w_stress
