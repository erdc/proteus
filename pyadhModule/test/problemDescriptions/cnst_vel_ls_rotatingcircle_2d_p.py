from pyadh import *
from pyadh.default_p import *
from math import *
## \page Tests Test Problems 
# \ref cnst_vel_ls_rotatingcircle_2d_p.py "Constant velocity LS (rotating circle)"

##\ingroup test
# \file cnst_vel_ls_rotatingcircle_2d_p.py
#
# @{
#  \brief Interface propagation using a level-set description and Hamilton-Jacobi formulation 
#  with a rotating velocity field.
#
#  The initial interface is described by the signed distance to a
#  circle.
#  
# \f{eqnarray*}
# \phi_t + s_n\|\nabla \phi\| &=& 0 \\ 
# \Omega &=& [0,1] \times [0,1] \\
# s_n &=& \vec u \cdot \vec n \\
# \vec n &=& \frac{\nabla \phi}{\|\nabla \phi\|} \\
#  u^{x} &=&  2\pi(y-1/2)\\
#  u^{y} &=&  2\pi(1/2-x)\\
# \phi^{ex}(x,y,t) &=& \sqrt{(x-x_c)^2 + (y-y_c)^2} - r
# \f}
# Here \f$r\f$ is the circle radius, \f$ x_c = \sin(2\pi t)/4 + 1/2\f$, and \f$\; y_c = \cos(2\pi t)/4 + 1/2 \f$.
# By default, \f$r=1/8\f$.
# Outflow boundaries are applied on \f$\partial \Omega\f$.
#
# \image html  save_cnst_vel_ls_rotatingcircle_2d_c0p1_exact.jpg "exact solution, T=1.0"
# \image latex save_cnst_vel_ls_rotatingcircle_2d_c0p1_exact.eps "exact solution, T=1.0"
# \image html  save_cnst_vel_ls_rotatingcircle_2d_c0p1_phi.jpg "phi"
# \image latex save_cnst_vel_ls_rotatingcircle_2d_c0p1_phi.eps "phi"
#
# 

nd = 2



class RotatingVelocityCircle:
    def __init__(self,radius=0.1):
        self.radius = radius
    def uOfXT(self,x,t):
        centerX = 0.25*sin(2*pi*t) + 0.5
        centerY = 0.25*cos(2*pi*t) + 0.5
        return sqrt((x[0]-centerX)**2 + (x[1]-centerY)**2)-self.radius

r0 = 1./8.

analyticalSolution = {0:RotatingVelocityCircle(r0)}

                                          

coefficients = RotatingVelocityLevelSet()

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

T = 1.0

## @}
