from pyadh import *
from pyadh.default_p import *
from math import *

##\page Tests Test Problems
#\ref cnst_vel_ls_circle_2d_p.py "Constant velocity LS (circle)"
#

##\ingroup test
#\file cnst_vel_ls_circle_2d_p.py
#
# @{
#
# 
# \brief Interface propagation using a level-set description and Hamilton-Jacobi formulation 
#  with a specified velocity field.
#
#  The initial interface is described by the signed distance to a
#  circle.
#  
# \f{eqnarray*}
# \phi_t + s_n\|\nabla \phi\| &=& 0 \\ 
# \Omega &=& [0,1] \times [0,1] \\
# s_n &=& \vec u \cdot \vec n \\
# \vec n &=& \frac{\nabla \phi}{\|\nabla \phi\|} \\
# \vec u &=& \vec u_0 \\
# \phi^{ex}(x,y,t) &=& \sqrt{(x-x_c)^2 + (y-y_c)^2} - r
# \f}
# Here \f$r\f$ is the circle radius, and \f$\vec x_c = \vec u t + \vec x_c^0 \f$.
# By default \f$r=1/8\f$, \f$\vec x_c^0 = (1/4,3/4)^T\f$, and \f$\vec u_0=(1,-1)^T\f$.
# Outflow boundaries are applied on \f$\partial \Omega\f$.
#
# \image html  save_cnst_vel_ls_circle_2d_c0p1_exact.jpg "exact solution, T=1.0"
# \image latex save_cnst_vel_ls_circle_2d_c0p1_exact.eps "exact solution, T=1.0"
# \image html  save_cnst_vel_ls_circle_2d_c0p1_phi.jpg "phi"
# \image latex save_cnst_vel_ls_circle_2d_c0p1_phi.eps "phi"
#
# 
#
#


nd = 2



class ConstantVelocityCircle:
    def __init__(self,radius=0.1,b=[1.,1.],startX=0.25,startY=0.5):
        self.radius = radius
        self.b      = b
        self.startX = startX
        self.startY = startY
    def uOfXT(self,x,t):
        centerX = self.b[0]*t + self.startX
        centerY = self.b[1]*t + self.startY
        return sqrt((x[0]-centerX)**2 + (x[1]-centerY)**2)-self.radius

x0 = 0.25
y0 = 0.75
r0 = 1./8.
b0 = Numeric.array([1.,-1.0],Numeric.Float)

analyticalSolution = {0:ConstantVelocityCircle(r0,b0,x0,y0)}
                                          
coefficients = ConstantVelocityLevelSet(b0)


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
#fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 0.5 #0.75

## @}
