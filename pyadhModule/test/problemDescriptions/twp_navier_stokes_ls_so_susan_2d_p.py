from pyadh import *
from pyadh.default_p import *
from susan import *

"""
Two-phase incompressible Navier-Stokes flow of waves around a ship hull
"""

## \page Tests Test Problems
#\ref twp_navier_stokes_ls_so_susan_2d_p.py "Two-phase incompressible Navier-Stokes flow of waves around a ship hull"
# 

##\ingroup test
# \file twp_navier_stokes_ls_so_susan_2d_p.py
# @{
#
# \brief Two-phase incompressible Navier-Stokes flow of waves around a ship hull
#
# The governing equations are described by the
# pyadh::TransportCoefficients::TwophaseNavierStokes_LS_SO class. The
# velocities are specified on he right hand boundary to reflect the
# ship's speed. The domain and initial/boundary conditions elsewhere
# are slip conditions and no flow on the lower boundary, constant
# pressure/outflow on the top boundary and outlfow on the left
# boundary.
#
#\image html susan_1.jpg
#\image html susan_2.jpg
#
# <A href="https://adh.usace.army.mil/pyadh-images/susan.avi"> AVI Animation of free surface and pressure </A>

                     
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

