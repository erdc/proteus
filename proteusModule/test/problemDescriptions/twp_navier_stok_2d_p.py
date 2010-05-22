from pyadh import *
from pyadh.default_p import *
from sloshbox import *

                     
analyticalSolution = None

#coefficients = TwophaseStokes_LS_SO(epsFact=1.5,g=[0.0,-9.8],nd=nd,steady=False)
#coefficients = TwophaseStokes_LS_SO(epsFact=0.5,g=[0.0,-9.8],nd=nd,steady=True)
coefficients = TwophaseNavierStokes_LS_SO(epsFact=epsFact,
                                          rho_0 = rho_0,
                                          nu_0 = nu_0,
                                          rho_1 = rho_1,
                                          nu_1 = nu_1,
                                          g=g,
                                          nd=nd)
#coefficients = TwophaseStokes_LS_SO(g=[0.0,9.8],nd=nd,steady=True)
#coefficients.rho_1 = coefficients.rho_0
#coefficients.mu_1 = coefficients.mu_0

