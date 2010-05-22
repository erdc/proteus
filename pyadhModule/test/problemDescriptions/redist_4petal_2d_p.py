from pyadh import *
from pyadh.default_p import *
from math import *

## \page Tests Test Problems 
# \ref redist_4petal_2d_p.py "Hamilton Jacobi LS redistancing (petal)"
# 

##\ingroup test
# \file redist_4petal_2d_p.py
# @{
#
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
# In this case, the fixed interface has four ''petals.'' The initial
# condition is the squared signed distance to the interface. The
# four-petal interface is given in polar coordinates as
#
# \f{eqnarray*}
# r(\theta) &=& r_0 + \frac{\cos(b\theta)}{a r_0} 
# \f}
# with \f$r_0 = 1/4\f$, \f$a = 40\f$, and \f$b = 4\f$, by default. 
#
# \image html  save_redist_4petal_2d_c0p1_phi0.jpg "initial condition"
# \image latex save_redist_4petal_2d_c0p1_phi0.eps "initial condition"
# \image html  save_redist_4petal_2d_c0p1_phi.jpg "signed distance solution"
# \image latex save_redist_4petal_2d_c0p1_phi.eps "signed distance solution"
#


nd = 2

#graph of curve in polar coords is
#r(\theta)=r_0 + cos(b\theta)/(a r_0)
#
#with r_0 = 0.25, a = 40, b=4
class ConstantVelocity4PetalSquared:
    def __init__(self,r0=0.15,b=4,a=10):
        self.r0     = r0
        self.b      = b
        self.a      = a
    def uOfXT(self,x,t):
        tx= x[0]-0.5; ty=x[1]-0.5
        r = sqrt(tx**2 + ty**2)
        th= atan2(tx,ty)
        pr= 0.5*(self.r0 + cos(self.b*th)/(self.a*self.r0))
        d = r-pr
        d2= r**2 - pr**2
        return d2
class ConstantVelocity4Petal:
    def __init__(self,r0=0.15,b=4,a=10):
        self.r0     = r0
        self.b      = b
        self.a      = a
    def uOfXT(self,x,t):
        tx= x[0]-0.5; ty=x[1]-0.5
        r = sqrt(tx**2 + ty**2)
        th= atan2(tx,ty)
        pr= 0.5*(self.r0 + cos(self.b*th)/(self.a*self.r0))
        d = r-pr
        d2= r**2 - pr**2
        return d
r0=0.25
a= 40
b=4
analyticalSolution = {0:ConstantVelocity4Petal(r0=r0,b=b,a=a)}

def Heaviside(p,eps=1.0e-2):
    if p < -eps:
        return 0.0
    if p > eps:
        return 1.0
    return 0.5*(1.+p/eps + 1./math.pi*math.sin(math.pi*p/eps))
def Se(p,eps=1.0e-2):
    return 2.0*(Heaviside(p,eps)-0.5)

                    

#coefficients = RedistanceLevelSetOld()
coefficients = RedistanceLevelSet(epsFact=1.5,u0=ConstantVelocity4PetalSquared())

#now define the Dirichlet boundary conditions

def getDBC(x,flag):
    pass
    
dirichletConditions = {0:getDBC}

weakDirichletConditions = {0:coefficients.setZeroLSweakDirichletBCs}#{0:getZeroLSDBC}
#weakDirichletConditions = None

initialConditions  = {0:ConstantVelocity4PetalSquared(r0=r0,b=b,a=a)}

fluxBoundaryConditions = {0:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 20.0e0

## @}
