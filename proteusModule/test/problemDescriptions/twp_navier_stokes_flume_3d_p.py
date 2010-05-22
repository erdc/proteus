from pyadh import *
from pyadh.default_p import *
from flume3d import *
"""
Two-phase incompressible Navier-Stokes flow of an air flume in water
"""

## \page Tests Test Problems
#\ref twp_navier_stokes_ls_so_flume_3d_p.py "Two-phase incompressible Navier-Stokes flow of an air flume in water"
# 

##\ingroup test
# \file twp_navier_stokes_ls_so_flume_3d_p.py
# @{
#
# \brief Two-phase incompressible Navier-Stokes flow of an air flume in water.
#
# The governing equations are described by the
# pyadh::TransportCoefficients::TwophaseNavierStokes_LS_SO class. The
# domain is the unit square and the boundary conditions are slip and
# no flow on the whole boundary. The pressure is set (arbitrarily) to
# zero at the node in the top left corner.
#
#
#\image html flume_1.jpg
#\image html flume_2.jpg
#
# <A href="https://adh.usace.army.mil/pyadh-images/flume.avi"> AVI Animation of velocity and pressure </A>
#

                     
analyticalSolution = None

coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=epsFact_viscosity,
                                             sigma=sigma_01,
                                             rho_0 = rho_0,
                                             nu_0 = nu_0,
                                             rho_1 = rho_1,
                                             nu_1 = nu_1,
                                             g=g,
                                             nd=nd,
                                             LS_model=2,
                                             KN_model=0,
                                             epsFact_density=epsFact_density,
                                             stokes=useStokes)
