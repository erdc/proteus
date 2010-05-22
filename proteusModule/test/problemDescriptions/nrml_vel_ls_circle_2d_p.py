from pyadh import *
from pyadh.default_p import *
from math import *

## \page Tests Test Problems 
# \ref nrml_vel_ls_circle_2d_p.py "Hamilton Jacobi LS normal speed propagation (circle)"
# 

##\ingroup test
# \file nrml_vel_ls_circle_2d_p.py
# @{
#  \brief Propagation of an interface in its normal direction using
#  a level-set description and Hamilton-Jacobi formulation.
#
#  The initial interface is described by the squared signed distance to a circle. 
#
# \f{eqnarray*}
# \phi_t + s_n\|\nabla \phi\| &=& 0 \\ 
# \Omega &=& [0,1] \times [0,1] \\
# s_n &=& 1 \\
# \vec n &=& \frac{\nabla \phi}{\|\nabla \phi\|} \\
# \phi^{0} &=& \sqrt{\left(x-x_c\right)^2 + \left(y-y_c\right)^2} - r 
# \f}
# By default, the circle center is \f$\vec x_c = (1/2,1/2)^T\f$, and the radius,\f$r=1/8\f$.
# Outflow boundaries are applied on \f$\partial \Omega\f$.
#
# \image html  save_nrml_vel_ls_circle_2d_c0p1_phi0.jpg "initial condition"
# \image latex save_nrml_vel_ls_circle_2d_c0p1_phi0.eps "initial condition"
# \image html  save_nrml_vel_ls_circle_2d_c0p1_phi.jpg "phi, T=0.25"
# \image latex save_nrml_vel_ls_circle_2d_c0p1_phi.eps "phi, T=0.25"


nd = 2



class UnitNormalVelocityCircle:
    def __init__(self,radius=0.1,startX=0.5,startY=0.5):
        self.radius = radius
        self.startX = startX
        self.startY = startY
    def uOfXT(self,x,t):
        centerX = self.startX
        centerY = self.startY
        return sqrt((x[0]-centerX)**2 + (x[1]-centerY)**2)-self.radius

analyticalSolution = {0:UnitNormalVelocityCircle(1./8.,0.5,0.5)}

class UnitNormalVelocityLevelSet(TransportCoefficients.TC_base):
    from pyadh.ctransportCoefficients import constantNormalVelocityLevelSetEvaluate
    def __init__(self,bn=1.0):
        self.bn=bn
        mass={0:{0:'linear'}}
#mwf need for grad(u)
        advection={0:{0:'linear'}}
        diffusion={}
        potential={}
#        reaction={0:{0:'constant'}}
        reaction={}
        hamiltonian={0:{0:'nonlinear'}}
        TransportCoefficients.TC_base.__init__(self,
                                             1,
                                             mass,
                                             advection,
                                             diffusion,
                                             potential,
                                             reaction,
                                             hamiltonian)
    def evaluate(self,t,c):
        self.constantNormalVelocityLevelSetEvaluate(self.bn,
                                                    c['x'],
                                                    c[('u',0)],
                                                    c[('grad(u)',0)],
                                                    c[('m',0)],
                                                    c[('dm',0,0)],
                                                    c[('f',0)],
                                                    c[('df',0,0)],
                                                    c[('H',0)],
                                                    c[('dH',0,0)]);


        #mwf when using r term in general HJ formulation
        if ('r',0) in c.keys():
            c[('r',0)].flat[:] = 0.0
#         if ('m',0) in c.keys():
#             c[('m',0)].flat[:] = c[('u',0)].flat[:]
#         if ('dm',0,0) in c.keys():
#             c[('dm',0,0)].flat[:] = 1.0
#         if ('f',0) in c.keys():
#             c[('f',0)].flat[:] = 0.0
#         if ('df',0,0) in c.keys():
#             c[('df',0,0)].flat[:] = 0.0
#         for i in range(len(c[('u',0)].flat)):
#             normgrad=sqrt(Numeric.dot(c[('grad(u)',0)].flat[nd*i:nd*(i+1)],
#                                       c[('grad(u)',0)].flat[nd*i:nd*(i+1)]))
#             c[('H',0)].flat[i]= normgrad
#             c[('dH',0,0)].flat[nd*i:nd*(i+1)]=c[('grad(u)',0)].flat[nd*i:nd*(i+1)]/normgrad
#         #end for
    #end def
                                          

coefficients = UnitNormalVelocityLevelSet()

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

T = 2.5e-1 #5.0e-1

## @}
