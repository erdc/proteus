from pyadh import *
from pyadh.default_p import *
from math import *

## \page Tests Test Problems 
# \ref redist_circle_2d_p.py "Hamilton Jacobi LS redistancing (circle)"
# 

##\ingroup test
# \file redist_circle_2d_p.py
# @{
#  \brief Solve the level set redistancing equation for a fixed interface. 
#
#  The redistancing function is the steady state solution of
# \f{eqnarray*}
# \phi_t + S\|\nabla \phi\| &=& S \\ 
# \Omega &=& [0,1] \times [0,1] \\
#  S &=& 2\left(H(\phi^0)-\frac{1}{2}\right) \\
# \f}
# where \f$H(u)\f$ is the Heaviside function. In the implementation, a smoothed
# version is used
# \f{eqnarray*}
# H_{\epsilon}(u) = \left\{ \begin{array}{ll}
# 0 & u < -\epsilon \\
# 1 & u > \epsilon  \\
# \frac{1}{2}(1 + \frac{u}{\epsilon} + \frac{1}{\pi}\sin(\frac{\pi u}{\epsilon})) & \mbox{otherwise}
# \end{array} \right.
# \f}
#
# In this case, the fixed interface is a circle. The initial condition is the
# squared signed distance to the circle. By default, the initial condition is
#
# \f{eqnarray*}
# \phi^0 &=& \left(x-\frac{1}{2}\right)^2 + \left(y - \frac{1}{2}\right)^2 - \frac{1}{16} 
# \f}
#
# \image html  save_redist_circle_2d_c0p1_phi0.jpg "initial condition"
# \image latex save_redist_circle_2d_c0p1_phi0.eps "initial condition"
# \image html  save_redist_circle_2d_c0p1_phi.jpg "signed distance solution"
# \image latex save_redist_circle_2d_c0p1_phi.eps "signed distance solution"
#
#
#

nd = 2

name = "circle_fmm_pwl_bdf2"
#mwf hack test unstructured mesh
#polyfile = "unstSquare"
from pyadh import RDLS
tryRDLS = True
if tryRDLS:
    LevelModelType = RDLS.OneLevelRDLS


class ConstantVelocityCircleSquared:
    def __init__(self,radius=0.1,b=[1.,1.],startX=0.25,startY=0.5):
        self.radius = radius
        self.b      = b
        self.startX = startX
        self.startY = startY
    def uOfXT(self,x,t):
        centerX = self.b[0]*0.0 + self.startX
        centerY = self.b[1]*0.0 + self.startY
        return (x[0]-centerX)**2 + (x[1]-centerY)**2-self.radius**2

class ConstantVelocityCircle:
    def __init__(self,radius=0.1,b=[1.,1.],startX=0.25,startY=0.5):
        self.radius = radius
        self.b      = b
        self.startX = startX
        self.startY = startY
    def uOfXT(self,x,t):
        centerX = self.b[0]*0.0 + self.startX
        centerY = self.b[1]*0.0 + self.startY
        return sqrt((x[0]-centerX)**2 + (x[1]-centerY)**2)-self.radius

x0 = 0.50
y0 = 0.50
r0 = 1./8.
b0 = numpy.array([-1.,0],'d')

analyticalSolution = {0:ConstantVelocityCircle(r0,b0,x0,y0)}


coefficients = RedistanceLevelSet(epsFact=0.25,u0=ConstantVelocityCircleSquared(r0,b0,x0,y0))

coefficients.variableNames=['u']
freezeLevelSet = False#True
if freezeLevelSet:
    if LevelModelType == RDLS.OneLevelRDLS:
        weakDirichletConditions = {0:RDLS.setZeroLSweakDirichletBCsSimple}
    else:
        weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}

#now define the Dirichlet boundary conditions

def getDBC(x,flag):
    pass
    #if (x[1] == 0.0):
    #    return lambda x,t: 0.0
    #if (x[0] == 0.0 or
    #    x[0] == 1.0 or
    #    x[1] == 0.0 or
    #    x[1] == 1.0):
    #    return lambda x,t: 0.0
    
dirichletConditions = {0:getDBC}

initialConditions  = {0:ConstantVelocityCircleSquared(r0,b0,x0,y0)}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 0.75e1

## @}
