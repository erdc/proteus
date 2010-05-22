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
#center_Y = 1.15
#center_X = 0.35
center_X = 0.5
center_Y = -0.5
radius = 1.0

eps=0.0#1.0e-1
#eps=0.01

def phi_circularLevelSetsDBC(x,t):
    import math
    import random
    random.seed(0.1)
    X_0 = x[0] - center_X
    Y_0 = x[1] - center_Y
    phi_0 = math.sqrt(X_0**2 + Y_0**2)
    return (1.0 + random.uniform(-eps,eps))*(phi_0 - radius)
    #return (1.0 + random.uniform(-eps,eps))*(x[1] - 0.5*L[1])
    #return x[0] - 0.5*L[0]
    #return x[0] - 0.5*L[0]+ x[1] - 0.5*L[1]
def phi_circularLevelSets(x,phi):
    import math
    import random
    nPoints=1
    for d in (x.shape[:-1]):
        nPoints*=d
    random.seed(0.1)
    for i in range(nPoints):
        X_0 = x.flat[i*3+0] - center_X
        Y_0 = x.flat[i*3+1] - center_Y
        phi_0 = math.sqrt(X_0**2 + Y_0**2)
        phi.flat[i]   = (1.0+random.uniform(-eps,eps))*(phi_0 - radius)
        #phi.flat[i] = (1.0+random.uniform(-eps,eps))*(x.flat[i*3+1] - 0.5*L[1])
#         phi.flat[i] = x.flat[i*3+0] - 0.5*L[0]
        #phi.flat[i] = x.flat[i*3+0] - 0.5*L[0] + x.flat[i*3+1] - 0.5*L[1]
#         if phi_0 < 1.0e-8:
#             phi_0 = 1.0e-8
#         grad_phi.flat[i*2 + 0] = X_0/phi_0
#         grad_phi.flat[i*2 + 1] = Y_0/phi_0
#         grad_phi.flat[i*2 + 0] += random.uniform(-eps,eps)
#         grad_phi.flat[i*2 + 1] += random.uniform(-eps,eps)


coefficients = LevelSetNormalCoefficients(epsFact=0.0,phi_func=phi_circularLevelSets)

#dummy to take out boundary nodes from problem
# def getDBC(x):
#    if (x[0] == 0.0 or
#        x[0] == L[0] or
#        x[1] == 0.0 or
#        x[1] == L[1]):
#        return phi_circularLevelSetsDBC

def getDBC(x):
     pass

dirichletConditions = {0:getDBC}

#fluxBoundaryConditions = {0}
## @}
