from pyadh import *
from pyadh.default_p import *
from math import *

## \page Tests Test Problems 
# \ref curvature_4petal_2d_p.py "Solve for the  level set curvature (petal problem)"
# 

##\ingroup test
# \file curvature_4petal_2d_p.py
# @{
#
#  \brief Solve the level set curvature equation for a fixed interface. 
#
#  The curvature function is the solution of
# \f{eqnarray*}
# \kappa = \nabla \cdot \nabla \phi \\
# \Omega &=& [0,1] \times [0,1] \\
# \f}
# where \f$\nabla \phi\f$ is the gradient of a given level set function
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

center_X = -0.15
center_Y = 0.35
center_Y = -0.15
center_X = 0.35
center_X = 1.15
center_Y = 0.35
center_Y = 0.5
center_X = 1.5
radius = 0.5

def grad_phi_circularLevelSets(x,grad_phi):
    import math
    import random
    grad_phi.flat[:]=-1
    nPoints=1
    for d in (x.shape[:-1]):
        nPoints*=d
    random.seed(0.1)
    eps=1.0e-2
    for i in range(nPoints):
        X_0 = x.flat[i*3+0] - center_X
        Y_0 = x.flat[i*3+1] - center_Y
        phi_0 = math.sqrt(X_0**2 + Y_0**2)
        phi   = phi_0 - radius
        if phi_0 < 1.0e-8:
            phi_0 = 1.0e-8
        grad_phi.flat[i*2 + 0] = X_0/phi_0
        grad_phi.flat[i*2 + 1] = Y_0/phi_0
#         grad_phi.flat[i*2 + 0] += random.uniform(-eps,eps)
#         grad_phi.flat[i*2 + 1] += random.uniform(-eps,eps)

class curvature_circularLevelSets:
    def uOfXT(self,x,t):
        import math
        X_0 = x[0] - center_X
        Y_0 = x[1] - center_Y
        phi_0 = math.sqrt(X_0**2 + Y_0**2)
        if phi_0 < 1.0e-8:
            phi_0 = 1.0e-8
        return -(2.0/phi_0 - X_0**2/phi_0**3 - Y_0**2/phi_0**3)

center_Xp = 0.5
center_Yp = 0.5
r0 = 0.5
b = 4
a = 10


theta=pi/4.0

n = [cos(theta),sin(theta)]

def grad_phi_straightLevelSets(x,grad_phi):
    import math
    import random
    grad_phi.flat[:]=-1
    nPoints=1
    for d in (x.shape[:-1]):
        nPoints*=d
    #random.seed(0.1)
    #eps=1.0e-3
    for i in range(nPoints):
        grad_phi.flat[i*2 + 0] = n[0]
        grad_phi.flat[i*2 + 1] = n[1]
        #grad_phi.flat[i*2 + 0] += random.uniform(-eps,eps)
        #grad_phi.flat[i*2 + 1] += random.uniform(-eps,eps)

class curvature_straightLevelSets:
    def uOfXT(self,x,t):
        return 0.0


analyticalSolution = {0:curvature_circularLevelSets()}
#analyticalSolution = {0:curvature_straightLevelSets()}
#analyticalSolution = None

coefficients = LevelSetCurvatureCoefficients()
coefficients = LevelSetCurvatureCoefficients(grad_phi_func=grad_phi_circularLevelSets)#grad_phi_func=grad_phi_4PetalSquared)
#coefficients = LevelSetCurvatureCoefficients(grad_phi_func=grad_phi_straightLevelSets)#grad_phi_func=grad_phi_4PetalSquared)
#coefficients = LevelSetCurvatureCoefficients(grad_phi_func=grad_phi_4petalLevelSets)#grad_phi_func=grad_phi_4PetalSquared)
fluxBoundaryConditions = {0:'outFlow'}
def getDBC(x):
    pass

dirichletConditions = {0:getDBC}
## @}
