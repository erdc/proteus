from pyadh import *
from pyadh.default_p import *
from math import *

## \page Tests Test Problems 
# \ref la_vortex_2d_p.py "Linear advection of a circular level set function in an oscillating vortex velocity field"
# 

##\ingroup test
# \file la_vortex_2d_p.py
# @{
#  \brief Conservative linear advection of a circle signed distance function
#  in a oscillating vortex velocity field.
#  
# \f{eqnarray*}
# \phi_t + \nabla \cdot (\vec u \phi) &=& 0 \\ 
# \Omega &=& [0,1] \times [0,1] \\
#  u^{x} &=& \cos(\pi t/8)\sin(2\pi y)\sin^2(\pi x) \\
#  u^{y} &=& -\cos(\pi t/8)\sin(2\pi x)\sin^{2}(\pi y) \\
#  \phi^{0} &=& \left(x-\frac{1}{2}\right)^2 + \left(y-\frac{3}{4}\right)^2 - 0.15^2
# \f}
# The solution should return to the initial condition at \f$T=8\f$.
# Outflow boundaries are applied on \f$\partial \Omega\f$.
# 
#
# \image html  save_la_vortex_2d_dgp2_exact.jpg "exact solution, T=8.0"
# \image latex save_la_vortex_2d_dgp2_exact.eps "exact solution, T=8.0"
# \image html  save_la_vortex_2d_dgp2_phi.jpg "RKDG P^2 solution, Cr=0.1, L^2 error= 7.84e-3"
# \image latex save_la_vortex_2d_dgp2_phi.eps "RKDG $P^2$ solution, Cr=0.1, $L^2$ error= 7.84e-3"
#
from pyadh import LADRellam
LevelModelType = LADRellam.OneLevelLADR
name = "la_vortex_ellam_2d_test"

nd = 2
nc = 1
#dispersivities
alpha_L = 0.0
alpha_T = 0.0
#molecular diffusion
d       = 0.0
#porosity
omega   = 1.0

class OscillatingVortex2D:
    def __init__(self):
        self.radius = 0.15
        self.xc=0.5
        self.yc=0.75
    def uOfXT(self,x,t):
        return (x[0]-self.xc)**2 + (x[1]-self.yc)**2 - self.radius**2
        
analyticalSolution = {0:OscillatingVortex2D()}

class DarcyVelocity:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        one8 = 1.0/8.0; xk = x[0]; yk=x[1]
        vx = cos(pi*one8*t)*sin(2.0*pi*yk)*sin(pi*xk)*sin(pi*xk)
        vy =-cos(pi*one8*t)*sin(2.0*pi*xk)*sin(pi*yk)*sin(pi*yk)
        return numpy.array([vx,vy])
darcyVelocityEval = {0:lambda x,t: DarcyVelocity().uOfXT(x,t)}
analyticalSolutionVelocity = {0:DarcyVelocity()}
velocitySpaceFlag = 'c0p1'#'bdm1'#'c0p1'
velocityIsTransient= True
coefficients = SubsurfaceTransportCoefficients.GroundwaterTransportCoefficientsELLAM(nc=nc,
                                                                                     omega=omega,
                                                                                     alpha_L=alpha_L,
                                                                                     alpha_T=alpha_T,
                                                                                     d = d,
                                                                                     nd = nd,
                                                                                     velocityFunctions = darcyVelocityEval,
                                                                                     velocitySpaceFlag = velocitySpaceFlag,
                                                                                     flowIsTransient=velocityIsTransient)

coefficients.variableNames=['u']

#now define the Dirichlet boundary conditions

def getDBC(x,tag):
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
def zeroadv(x):
    return lambda x,t: 0.0

def getAFBC_ci(x,tag,ci=0):
    if x[0] == 0.0 and analyticalSolutionVelocity[ci].uOfXT(x,0.)[0] > 0.:
        return lambda x,t: 0.0
    if x[0] == 1.0 and analyticalSolutionVelocity[ci].uOfXT(x,0.)[0] < 0.:
        return lambda x,t: 0.0
    if x[1] == 0.0 and analyticalSolutionVelocity[ci].uOfXT(x,0.)[1] > 0.:
        return lambda x,t: 0.0
    if x[1] == 1.0 and analyticalSolutionVelocity[ci].uOfXT(x,0.)[1] < 0.:
        return lambda x,t: 0.0
advectiveFluxBoundaryConditions =  dict((cj,lambda x,tag,cii=cj : getAFBC_ci(x,tag,ci=cii)) for cj in range(nc))

def getDFBC_ci(x,tag,ci=0):
   if x[0] == 1.0 and analyticalSolutionVelocity[ci].uOfXT(x,0.)[0] > 0.:
       return lambda x,t: 0.0
   if x[0] == 0.0 and analyticalSolutionVelocity[ci].uOfXT(x,0.)[0] < 0.:
       return lambda x,t: 0.0
   if x[1] == 1.0 and analyticalSolutionVelocity[ci].uOfXT(x,0.)[1] > 0.:
       return lambda x,t: 0.0
   if x[1] == 0.0 and analyticalSolutionVelocity[ci].uOfXT(x,0.)[1] < 0.:
       return lambda x,t: 0.0



diffusiveFluxBoundaryConditions = {0:{}}
for ci in range(nc):
    diffusiveFluxBoundaryConditions[ci] = dict((ck,lambda x,tag,cii=ck : getDFBC_ci(x,tag,ci=cii)) for ck in range(nc))

T = 8.0

## @}
