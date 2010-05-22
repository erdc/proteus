from pyadh import *
from pyadh.default_p import *

## \page Tests Test Problems 
# \ref twp_darcy_ls_so_bubble_2d_p.py "Two-phase Darcy flow with a sharp interface approximation (bubble)"
# 

##\ingroup test
# \file twp_darcy_ls_so_bubble_2d_p.py
# @{
#
#  \brief Two-phase, incompressible flow using a sharp interface
#  approximation and assuming Darcy's law for momentum.
#
#  The interface is described using a level-set approximation.
#
# \f{eqnarray*}
# \nabla \cdot \vec q &=& 0 \\
# \vec q &=& -K\left(\nabla p - \rho\vec g\right) \\ 
# \rho &=& \rho_0 + H(\phi)\left(\rho_1 - \rho_0\right) \\
# K    &=& K_0 + H(\phi)\left(K_1 - K_0\right) 
# \f}
#
# Here, \f$H(u)\f$ is the Heaviside function, \f$\phi\f$ is the level
# set function defining the fluid/fluid interface, and the jump
# conditions \f$[p] = 0\f$ and \f$[\vec q\cdot \vec n]=0\f$ are used
# implicitly.
#
# In the implementation, a smoothed version of the Heaviside function is used
# \f{eqnarray*}
# H_{\epsilon}(u) = \left\{ \begin{array}{ll}
# 0 & u < -\epsilon \\
# 1 & u > \epsilon  \\
# \frac{1}{2}(1 + \frac{u}{\epsilon} + \frac{1}{\pi}\sin(\frac{\pi u}{\epsilon})) & \mbox{otherwise}
# \end{array} \right.
# \f}
#
# The level set equation for interface propagation used is
# \f{eqnarray*}
# \phi_t + \vec q \cdot \nabla \phi &=& 0 \\ 
# \f}
#
# The two fluids are air and water, so that \f$\rho_0/\rho_1 = 10^{-3}\f$ and \f$K_0/K_1 = 55\f$.
# By default, the domain is \f$\Omega = [0,1] \times [0,1] \f$ and the boundary conditions
# are hydrostatic for a given phase,
# while the initial condition is a circular bubble of radius \f$1/8\f$ 
#
# \f{eqnarray*}
# \phi^0(x,y) &=& \sqrt{\left(x-x_c\right)^2 + \left(y-y_c\right)^2} - \frac{1}{8} \\
#  p^0(x,y) &=& \rho_1(1-y), \mbox{ or }\\
# \phi^0(x,y) &=& -\sqrt{\left(x-x_c\right)^2 + \left(y-y_c\right)^2} + \frac{1}{8} \\
#  p^0(x,y) &=& \rho_0(1-y),\\
# \f}
#
# 
#
# \image html save_twp_darcy_ls_so_bubble_2d_c0p1_phi.jpg "level set solution"
# \image latex save_twp_darcy_ls_so_bubble_2d_c0p1_phi.eps "level set solution"
# \image html save_twp_darcy_ls_so_bubble_2d_c0p1_psi.jpg "p"
# \image latex save_twp_darcy_ls_so_bubble_2d_c0p1_psi.eps "p"
# \image html save_twp_darcy_ls_so_bubble_2d_c0p1_vel.jpg "velocity "
# \image latex save_twp_darcy_ls_so_bubble_2d_c0p1_vel.eps "velocity "
#
#

nd = 2
Ident = Numeric.zeros((nd,nd),Numeric.Float)
Ident[0,0]=1.0; Ident[1,1] = 1.0
gravity = Numeric.zeros((nd,),Numeric.Float)
gravity[1] = -1.0

rhoP = 1.0; rhoM = 1.e-3; #1.e-3
Kp   = 1.0; Km   = 5.5e1; #55.

heaviEps = 1.0e-2 #1e-2 worked ok w/o redist

class DarcyFlowLS(TransportCoefficients.TC_base):
    from pyadh.ctransportCoefficients import darcySharpInterfaceFlowEvaluate
    from pyadh.ctransportCoefficients import ncLevelSetCoefficientsEvaluate
    def __init__(self,Km,Kp,rhoM,rhoP,gravity,eps,LS_model=1):
        nc = 1
        self.Km = Km
        self.Kp = Kp
        self.rhoM = rhoM
        self.rhoP = rhoP
        self.gravity = gravity
        self.heaviEps = eps
        self.LS_model = LS_model
        mass     = {}
        advection = {}
        diffusion = {}
        potential = {}
        reaction  = {}
        hamiltonian = {}
        #Darcy mass balance
        i = 0
        mass[i] = {i : 'linear'}
        diffusion[i] = {i : {i:'constant'}}
        potential[i] = {i : 'u'}
        advection[i] = {i : 'constant'}
        reaction[i] = {i : 'constant'}
        TransportCoefficients.TC_base.__init__(self,
                                               nc,
                                               mass,
                                               advection,
                                               diffusion,
                                               potential,
                                               reaction,
                                               hamiltonian)
        self.useC = True
        self.q_uls  = None
        self.ebq_uls= None
    #end init
    def attachModels(self,modelList):
        self.q_uls = modelList[self.LS_model].q[('u',0)]
        self.ebq_uls = modelList[self.LS_model].ebq[('u',0)]
    def initializeElementQuadrature(self,t,cq):
        #initialize so it can run without a flow model
        self.q_uls = Numeric.ones(cq[('u',0)].shape,Numeric.Float)
    def initializeElementBoundaryQuadrature(self,t,cebq,cebq_global):
        self.ebq_uls = Numeric.ones(cebq[('u',0)].shape,Numeric.Float)
    
    def evaluate(self,t,c):
        if self.q_uls.shape == c[('u',0)].shape:
            uls = self.q_uls
        else:
            uls = self.ebq_uls
        self.darcySharpInterfaceFlowEvaluate(self.Km,self.rhoM,
                                             self.Kp,self.rhoP,
                                             self.heaviEps,
                                             self.gravity,
                                             c[('u',0)],
                                             c[('grad(u)',0)],
                                             uls,#level set function
                                             c[('phi',0)],#potential
                                             c[('a',0,0)],
                                             c[('f',0)],
                                             c[('r',0)],
                                             c[('m',0)],
                                             c[('dphi',0,0)],
                                             c[('da',0,0,0)],
                                             c[('df',0,0)],
                                             c[('dr',0,0)],
                                             c[('dm',0,0)])
    #eval
                
class HydrostaticSinglePhase_psi:
    def __init__(self,rho=1.0):
        self.rho  =rho
    def uOfX(self,X):
        uHydro = 0.0 + self.rho*(1.0+gravity[1]*X[1])
        return uHydro
    def uOfXT(self,X,t):
        return self.uOfX(X)
    #end
#end HydrostaticSinglePhase_psi


def getDBCHydroP(X):
    #background P fluid
    if X[0] in [0.0,1.0] or X[1] in [0.0,1.0]:
        return lambda x,t: HydrostaticSinglePhase_psi(rhoP).uOfXT(x,t)
#dirichlet
def getDBCHydroM(X):
    #background M fluid
    if X[0] in [0.0,1.0]:
        return lambda x,t: HydrostaticSinglePhase_psi(rhoM).uOfXT(x,t)
#dirichlet

analyticalSolution = None
analyticalSolutionVelocity = None

#background P fluid
#initialConditions  = {0:HydrostaticSinglePhase_psi(rhoP)}
#dirichletConditions= {0:getDBCHydroP}
#background M fluid
initialConditions  = {0:HydrostaticSinglePhase_psi(rhoM)}
dirichletConditions= {0:getDBCHydroM}

coefficients = DarcyFlowLS(Km,Kp,rhoM,rhoP,gravity,heaviEps)

coefficients.variableNames=['psi']
fluxBoundaryConditions = {0:'noFlow'}
advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}


## @}

