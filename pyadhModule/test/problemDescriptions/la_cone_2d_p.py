from pyadh import *
from pyadh.default_p import *
from math import *
"""
The 2D rotating gaussian cone problem docstring
"""

## \page Tests Test Problems 
# \ref la_cone_2d_p.py "Linear Advection (cone)"
#  

##\ingroup test
#  \file la_cone_2d_p.py
# @{
#
#  \brief Conservative linear advection of a cone in a rotating
#  velocity field.
#
# \f{eqnarray*}
# \phi_t + \nabla \cdot (\vec u \phi) &=& 0 \\ 
# \Omega &=& [0,1] \times [0,1] \\
#  u^{x} &=&  2\pi(y-1/2)\\
#  u^{y} &=&  2\pi(1/2-x)\\
#  \phi^{ex}(x,y,t) &=&  \left\{
#  \begin{array}{ll}
#  \left(1+\cos(\pi\bar{x}/r_0)\right)
#  \left(1+\cos(\pi\bar{y}/r_0)\right)/4, & \|\bar{\vec r}\| < r_0 \\
#  0 & \mbox{otherwise} \\
# \end{array} \right\}
#  \f}
#
# where \f$\bar{x} = x - x_c\f$, \f$\bar{y} = y - y_c\f$, and
# \f$ x_c = \sin(2\pi t)/4 + 1/2\f$, \f$\; y_c = \cos(2\pi t)/4 + 1/2 \f$
#
# \image html  save_la_cone_2d_dgp2_exact.jpg "exact solution, T=1.0"
# \image latex save_la_cone_2d_dgp2_exact.eps "exact solution, T=1.0"
# \image html  save_la_cone_2d_dgp2_phi.jpg "RKDG P^2 solution, Cr=0.15, L^2 error= 1.73e-3"
# \image latex save_la_cone_2d_dgp2_phi.eps "RKDG $P^2$ solution, Cr=0.15, $L^2$ error= 1.73e-3"
#

#
#
nd = 2

class RotatingCone2D:
    def __init__(self,radius):
        self.radius = radius
    def uOfXT(self,x,t):
        centerX = 0.25*sin(2*pi*t) + 0.5
        centerY = 0.25*cos(2*pi*t) + 0.5
        coneX = x[0] - centerX
        coneY = x[1] - centerY
        if sqrt(coneX**2 + coneY**2) < self.radius:
            return 0.25*(1.0 + cos(pi*coneX/self.radius))*(1.0+cos(pi*coneY/self.radius))
        else:
            return 0.0

analyticalSolution = {0:RotatingCone2D(1.0/8.0)}

coefficients = UnitSquareRotation()

coefficients.variableNames=['u']

#now define the Dirichlet boundary conditions

def getDBC(x):
    pass
    #dgp2 is ok with dirichlet bc's has problems on outflow?
    #if (x[0] == 0.0 or
    #    x[0] == 1.0 or
    #    x[1] == 0.0 or
    #    x[1] == 1.0):
    #    return lambda x,t: 0.0
#Recall
# v_x &=& 2 \pi (y - 1/2) \\
# v_y &=& 2 \pi (1/2-x) \\
def zeroInflow(x):
    if (x[0] == 0.0 and x[1] >= 0.5):
        return lambda x,t: 0.0
    if (x[0] == 1.0 and x[1] <= 0.5):
        return lambda x,t: 0.0
    if (x[1] == 0.0 and x[0] <= 0.5):
        return lambda x,t: 0.0
    if (x[1] == 1.0 and x[0] >= 0.5):
        return lambda x,t: 0.0
    
dirichletConditions = {0:getDBC}

initialConditions  = {0:analyticalSolution[0]}

fluxBoundaryConditions = {0:'outFlow'}

def zeroadv(x):
    return lambda x,t: 0.0
#advectiveFluxBoundaryConditions =  {}
advectiveFluxBoundaryConditions =  {0:zeroadv}

diffusiveFluxBoundaryConditions = {0:{}}

#mwf hack
T = 1.0

## @}
