from pyadh import *
from pyadh.default_p import *
from bubble3d import *
"""
Two-phase incompressible Navier-Stokes flow of an air bubble in water
"""

## \page Tests Test Problems
#\ref twp_navier_stokes_ls_so_bubble_2d_p.py "Two-phase incompressible Navier-Stokes flow of an air bubble in water"
# 

##\ingroup test
# \file twp_navier_stokes_ls_so_bubble_2d_p.py
# @{
#
# \brief Two-phase incompressible Navier-Stokes flow of an air bubble in water.
#
# The governing equations are described by the
# pyadh::TransportCoefficients::TwophaseNavierStokes_LS_SO class. The
# domain is the unit square and the boundary conditions are slip and
# no flow on the whole boundary. The pressure is set (arbitrarily) to
# zero at the node in the top left corner.
#
#
#\image html bubble_1.jpg
#\image html bubble_2.jpg
#
# <A href="https://adh.usace.army.mil/pyadh-images/bubble.avi"> AVI Animation of velocity and pressure </A>
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

dirichletConditions = {0:getDBC_p_bubble,
                       1:getDBC_u_bubble,
                       2:getDBC_v_bubble,
                       3:getDBC_w_bubble}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_bubble}

diffusiveFluxBoundaryConditions = {0:{}}

fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow',
                          3:'outFlow'}

initialConditions = {0:Hydrostatic_p(),
                     1:AtRest(),
                     2:AtRest(),
                     3:AtRest()}
