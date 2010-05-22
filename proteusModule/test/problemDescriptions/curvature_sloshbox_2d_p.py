from pyadh import *
from pyadh.default_p import *
from sloshbox import *
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

#coefficients = LevelSetCurvatureCoefficients(epsFact=0.01,LSModel_index=2)
coefficients = LevelSetCurvatureCoefficients(epsFact=epsFact_curvature,LSModel_index=2,nd=nd)

fluxBoundaryConditions = {0:'outFlow'}

def getDBC_kappa(x,flag):
    pass

dirichletConditions = {0:getDBC_kappa}

def getAFBC_kappa(x,flag):
     pass

advectiveFluxBoundaryConditions =  {0:getAFBC_kappa}

diffusiveFluxBoundaryConditions = {0:{}}
## @}
