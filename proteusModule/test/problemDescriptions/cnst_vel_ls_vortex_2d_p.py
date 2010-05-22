from pyadh import *
from pyadh.default_p import *
from math import *
## \page Tests Test Problems
# \ref cnst_vel_ls_vortex_2d_p.py "Constant velocity LS (vortex)"
#  

##\ingroup test
# \file cnst_vel_ls_vortex_2d_p.py
#
# @{
#  \brief Interface propagation using a level-set description and Hamilton-Jacobi formulation 
#  with a oscillating vortex velocity field.
#
#  The initial interface is described by the squared signed distance
#  to a circle.
#
# \f{eqnarray*}
# \phi_t + s_n\|\nabla \phi\| &=& 0 \\ 
# \Omega &=& [0,1] \times [0,1] \\
# s_n &=& \vec u \cdot \vec n \\
# \vec n &=& \frac{\nabla \phi}{\|\nabla \phi\|} \\
#  u^{x} &=& \cos(\pi t/8)\sin(2\pi y)\sin^2(\pi x) \\
#  u^{y} &=& -\cos(\pi t/8)\sin(2\pi x)\sin^{2}(\pi y) \\
# \phi^{0} &=& \left(x-\frac{1}{2}\right)^2 + \left(y-\frac{3}{4}\right)^2 - r^2
# \f}
# The circle radius,\f$r\f$, is \f$1/8\f$ by default.
# The solution should return to the initial condition at \f$T=8\f$.
# Outflow boundaries are applied on \f$\partial \Omega\f$.
#
# \image html  save_cnst_vel_ls_vortex_2d_c0p1_phi0.jpg "initial condition"
# \image latex save_cnst_vel_ls_vortex_2d_c0p1_phi0.eps "initial condition"
# \image html  save_cnst_vel_ls_vortex_2d_c0p1_phiT1_6.jpg "phi, T=1.6"
# \image latex save_cnst_vel_ls_vortex_2d_c0p1_phiT1_6.eps "phi, T=1.6"
#


nd = 2


class OscillatingVortex2D:
    def __init__(self):
        self.radius = 0.15
        self.xc=0.5
        self.yc=0.75
    def uOfXT(self,x,t):
        return (x[0]-self.xc)**2 + (x[1]-self.yc)**2 - self.radius**2


analyticalSolution = {0:OscillatingVortex2D()}


coefficients = UnitSquareVortexLevelSet()

coefficients.variableNames=['u']

#now define the Dirichlet boundary conditions

def getDBC(x):
    pass
    #if (x[1] == 0.0):
    #    return lambda x,t: 0.0
    #if (x[0] == 0.0 or
    #    x[0] == 1.0 or
    #    x[1] == 0.0 or
    #    x[1] == 1.0):
    #    return lambda x,t: 0.0
    
dirichletConditions = {0:getDBC}

initialConditions  = {0:analyticalSolution[0]}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 8.e0

## @}
