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
T = ct.T  # might not be necessary

LevelModelType = RANS2P.LevelModel
if ct.useOnlyVF:
    LS_model = None
else:
    LS_model = 2
if ct.useRANS >= 1:
    Closure_0_model = 5
    Closure_1_model = 6
    if ct.useOnlyVF:
        Closure_0_model = 2
        Closure_1_model = 3
    if ct.movingDomain:
        Closure_0_model += 1
        Closure_1_model += 1
else:
    Closure_0_model = None
    Closure_1_model = None

# for absorption zones (defined as regions)
if hasattr(domain, 'porosityTypes'):
    porosityTypes = domain.porosityTypes
    dragAlphaTypes = domain.dragAlphaTypes
    dragBetaTypes = domain.dragBetaTypes
    epsFact_solid = domain.epsFact_solid
else:
    porosityTypes = None
    dragAlphaTypes = None
    dragBetaTypes = None
    epsFact_solid = None

coefficients = RANS2P.Coefficients(epsFact=ct.epsFact_viscosity,
                                   sigma=0.0,
                                   rho_0=ct.rho_0,
                                   nu_0=ct.nu_0,
                                   rho_1=ct.rho_1,
                                   nu_1=ct.nu_1,
                                   g=ct.g,
                                   nd=nd,
                                   ME_model=int(ct.movingDomain)+0,
                                   VF_model=int(ct.movingDomain)+1,
                                   LS_model=int(ct.movingDomain)+LS_model,
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
                                   movingDomain=ct.movingDomain,
                                   porosityTypes=porosityTypes,
                                   dragAlphaTypes=dragAlphaTypes,
                                   dragBetaTypes=dragBetaTypes,
                                   epsFact_solid=epsFact_solid,
                                   barycenters=ct.domain.barycenters)


dirichletConditions = {0: lambda x, flag: domain.bc[flag].p_dirichlet.init_cython(),
                       1: lambda x, flag: domain.bc[flag].u_dirichlet.init_cython(),
                       2: lambda x, flag: domain.bc[flag].v_dirichlet.init_cython()}

advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].p_advective.init_cython(),
                                   1: lambda x, flag: domain.bc[flag].u_advective.init_cython(),
                                   2: lambda x, flag: domain.bc[flag].v_advective.init_cython()}

diffusiveFluxBoundaryConditions = {0: {},
                                   1: {1: lambda x, flag: domain.bc[flag].u_diffusive.init_cython()},
                                   2: {2: lambda x, flag: domain.bc[flag].v_diffusive.init_cython()}}

if nd == 3:
    dirichletConditions[3] = lambda x, flag: domain.bc[flag].w_dirichlet.init_cython()
    advectiveFluxBoundaryConditions[3] = lambda x, flag: domain.bc[flag].w_advective.init_cython()
    diffusiveFluxBoundaryConditions[3] = {3: lambda x, flag: domain.bc[flag].w_diffusive.init_cython()}

class P_IC:
    def uOfXT(self, x, t):
        return ct.twpflowPressure_init(x, t)

class U_IC:
    def uOfXT(self, x, t):
        return 0.0

class V_IC:
    def uOfXT(self, x, t):
        return 0.0

class W_IC:
    def uOfXT(self, x, t):
        return 0.0

initialConditions = {0: P_IC(),
                     1: U_IC(),
                     2: V_IC()}
if nd == 3:
    initialConditions[3] = W_IC()
