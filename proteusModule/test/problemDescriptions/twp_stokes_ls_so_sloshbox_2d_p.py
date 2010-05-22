from pyadh import *
from pyadh.default_p import *
from sloshbox import *
                     
analyticalSolution = None

coefficients = TwophaseStokes_LS_SO(epsFact = epsFact,
                                    rho_0 = rho_0,
                                    nu_0 = nu_0,
                                    rho_1 = rho_1,
                                    nu_1 = nu_1,
                                    g=g,
                                    nd=nd,
                                    steady=False)

