from pyadh import *
from pyadh.default_p import *
from wavetank import *
                     
coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=epsFact_viscosity,
                                             sigma=sigma_01,
                                             rho_0 = rho_0,
                                             nu_0 = nu_0,
                                             rho_1 = rho_1,
                                             nu_1 = nu_1,
                                             g=g,
                                             nd=nd,
                                             LS_model=1,
                                             KN_model=None,
                                             epsFact_density=epsFact_density,
                                             stokes=useStokes)

dirichletConditions = {0:getDBC_p_wavetank,
                       1:getDBC_u_wavetank,
                       2:getDBC_v_wavetank}

fluxBoundaryConditions = {0:'mixedFlow',
                          1:'mixedFlow',
                          2:'mixedFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_wavetank,
                                    1:getAFBC_u_wavetank,
                                    2:getAFBC_v_wavetank}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u_wavetank},
                                   2:{2:getDFBC_v_wavetank}}

initialConditions = {0:Hydrostatic_p(),
                     1:AtRest(),
                     2:AtRest()}
