from proteus.default_p import *
from proteus.mprans import RANS2P
import numpy as np
from proteus import Context

ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions

genMesh = mesh.genMesh
movingDomain = ct.movingDomain
T = ct.T  # might not be necessaryd

Closure_0_model = None
Closure_1_model = None

LevelModelType = RANS2P.LevelModel
coefficients = RANS2P.Coefficients(epsFact=ct.epsFact_viscosity,
                                   sigma=0.0,
                                   rho_0=ct.rho_0,
                                   nu_0=ct.nu_0,
                                   rho_1=ct.rho_1,
                                   nu_1=ct.nu_1,
                                   g=ct.g,
                                   nd=nd,
                                   ME_model=0,
                                   CLSVOF_model=ct.CLSVOF_model,
                                   VF_model=ct.VF_model,
                                   LS_model=ct.LS_model,
                                   Closure_0_model=Closure_0_model,
                                   Closure_1_model=Closure_1_model,
                                   epsFact_density=ct.epsFact_density,
                                   stokes=False,
                                   useVF=ct.useVF,
                                   useRBLES=ct.useRBLES,
                                   useMetrics=ct.useMetrics,
                                   eb_adjoint_sigma=1.0,
                                   eb_penalty_constant=ct.weak_bc_penalty_constant,
                                   forceStrongDirichlet=ct.ns_forceStrongDirichlet,
                                   turbulenceClosureModel=ct.ns_closure,
                                   movingDomain=ct.movingDomain)

    
dirichletConditions = {0: lambda x, flag: domain.bc[flag].p_dirichlet.init_cython(),
                       1: lambda x, flag: domain.bc[flag].u_dirichlet.init_cython(),
                       2: lambda x, flag: domain.bc[flag].v_dirichlet.init_cython()}

advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].p_advective.init_cython(),
                                   1: lambda x, flag: domain.bc[flag].u_advective.init_cython(),
                                   2: lambda x, flag: domain.bc[flag].v_advective.init_cython()}

diffusiveFluxBoundaryConditions = {0:{},
                                   1: {1: lambda x, flag: domain.bc[flag].u_diffusive.init_cython()},
                                   2: {2: lambda x, flag: domain.bc[flag].v_diffusive.init_cython()}}

class AtRest:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:AtRest(),
                     1:AtRest(),
                     2:AtRest()}
