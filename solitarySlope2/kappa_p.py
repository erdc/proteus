from proteus.default_p import *
from proteus.mprans import Kappa
from proteus import Context
ct = Context.get()
domain = ct.domain
nd = domain.nd

LevelModelType = Kappa.LevelModel
if ct.useOnlyVF:
    RD_model = None
    LS_model = None
    dissipation_model = 3
    ME_model = 2
else:
    RD_model = 3
    LS_model = 2
    ME_model = 5
    dissipation_model = 6
if ct.movingDomain:
    dissipation_model += 1
    ME_model += 1
    LS_model += 1
    RD_model += 1
#
dissipation_model_flag = 1
if ct.useRANS == 2:
    dissipation_model_flag = 2
elif ct.useRANS == 3:
    dissipation_model_flag = 3

coefficients = Kappa.Coefficients(V_model=int(ct.movingDomain)+0,
                                  ME_model=ME_model,
                                  LS_model=LS_model,
                                  RD_model=RD_model,
                                  dissipation_model=dissipation_model,
                                  dissipation_model_flag=dissipation_model_flag,#1 -- K-epsilon, 2 -- K-omega 1998, 3 -- K-omega 1988
                                  useMetrics=ct.useMetrics,
                                  rho_0=ct.rho_0,
                                  nu_0=ct.nu_0,
                                  rho_1=ct.rho_1,
                                  nu_1=ct.nu_1,
                                  g=ct.g,
                                  nd=nd,
                                  c_mu=0.09,
                                  sigma_k=1.0,
                                  sc_uref=ct.kappa_sc_uref,
                                  sc_beta=ct.kappa_sc_beta)


dirichletConditions = {0: lambda x, flag: domain.bc[flag].k_dirichlet.init_cython()}
advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].k_advective.init_cython()}
diffusiveFluxBoundaryConditions = {0: {0: lambda x, flag: domain.bc[flag].k_diffusive.init_cython()}}



class ConstantIC:
    def __init__(self, cval=0.):
        self.cval = cval
    def uOfXT(self, x, t):
        return self.cval

initialConditions = {0: ConstantIC(cval=0)}
