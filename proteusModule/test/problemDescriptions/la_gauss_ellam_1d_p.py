from pyadh import *
from pyadh.default_p import *
"""
1D,Linear advection of a Gaussian with ELLAM
"""

## \page Tests Test Problems 
# \ref la_gauss_ellam_1d_p.py "Linear advection-diffusion of a Gaussian with ELLAM"
#

##\ingroup test
#\file la_gauss_ellam_1d_p.py
#@{
#
# \brief Uncoupled, multicomponent linear advection-diffusion for gaussian initial conditions.
#
# The linear advection dispersion (reaction) equation is written here as
#
#\f{eqnarray*}
# u_t + v u_x -D u_xx - f &=& 0 \\
# v &=& q/omega \\
# D &=& \alpha_L|v| + d_m \\
# u(x,0) &=& e^{-(x-x_{c})^{2}/(2\sigma^{2})} \\
# u(x,t) &=& A e^{-B^2/C^2} \\
# A &=& \frac{2^{1/2}\sigma}{(2\sigma^2 + 4D*t)^{1/2}} \\
#   &=& \frac{2^{1/2}\sigma}{C} \\
# B &=& x - v*t - x_{c}
# f &=&  -2^{3/2}\sigma*C^{-3}*D*exp(-B^{2}/C^{2}) + 
#              A*[2B*v/C^{2} + 4*D*B^{2}/C^{4}]exp(-B^{2}/C^{2})
#      -v*A*2*B/C^{2}*exp( -B^{2}/C^{2})
#      -D*A*[-2/C^{2} + 4*B^{2}/C^{4}]*exp(-B^{2}/C^{2})
#\f}
#
# Here, $x_{c}$ is the center of the hump in the initial condition and 
#  $\sigma$ is the spread factor
# That the boundary conditions are 
#
#\f[
#  u(x_I,t) = 0.0,  -D\pd{u}{x} \cdot n = 0, x= x_O
# \f]
#where $x_I$ and $x_O$ are the inflow and outflow respectively
# Note that the coefficient class evaluates the equation in conservative form
#
# Below, component ci's velocity is scaled by (-1**(ci)) and the
# initial/boundary conditions are shifted accordingly to test uncoupled
# multicomponent implementation
# \todo 

from pyadh import LADRellam
LevelModelType = LADRellam.OneLevelLADR
name = "la_gauss_ellam_1d"

nd = 1
nc = 1
#dispersivities
alpha_L = 1.0e-4
alpha_T = 1.0e-5##immaterial in 1d
#molecular diffusion
d       = 0.0
#porosity
omega   = 0.5 #todo change so that input velocity is Darcy velocity not pore velocity
#constant and possible derivative in space for velocities
darcyVelocity_0 =  1.0
darcyVelocity_x =  0.5
velocity_0      = darcyVelocity_0/omega
velocity_x      = darcyVelocity_x/omega
velocitySpaceFlag = 'rt0'#'bdm1'#'c0p1'
velocityIsTransient=True
tForVelocityReversal = 0.25#0.25#10000.0
class DarcyVelocity:
    def __init__(self,v=darcyVelocity_0,vx=darcyVelocity_x,transient=velocityIsTransient,tForReversal=tForVelocityReversal):
        self.v = v; self.vx = vx
        self.transient=transient; self.tForReversal=tForReversal
    def uOfXT(self,x,t):
        if self.transient:
            v = (self.v + self.vx*x[0])*(self.tForReversal-t)/(self.tForReversal-0.0)
            return numpy.array([v])
        return numpy.array([self.v + self.vx*x[0]])
class ParticleVelocity:#aka adjointVelocity, for tracking tests
    def __init__(self,v=velocity_0,vx=velocity_x,transient=velocityIsTransient,tForReversal=tForVelocityReversal):
        self.v = v; self.vx = vx
        self.transient=transient; self.tForReversal=tForReversal
    def uOfXT(self,x,t):
        if self.transient:
            v = (self.v + self.vx*x[0])*(self.tForReversal-t)/(self.tForReversal-0.0)
            return numpy.array([v])
        return numpy.array([self.v + self.vx*x[0]])
class ParticleVelocityTrajectory:
    def __init__(self,v=velocity_0,vx=velocity_x):
        self.v = v; self.vx = vx
    def uOfXT(self,x,dt):
        if abs(self.vx) < 1.0e-10:
            return x[0] + dt*self.v
        else:
            v0 = self.v + self.vx*x[0]
            return x[0] + v0*(exp(self.vx*dt)-1.0)/self.vx
    def dtOfdX(self,x,dx,n=numpy.array([1.0])):
        v0 = n[0]*(self.v + self.vx*x[0])
        if abs(v0) < 1.0e-10:
            return 0.0
        if abs(self.vx) < 1.0e-10:
            return dx[0]/self.v
        if dx[0]*self.vx/v0 <= -1: return None
        return math.log(dx[0]*self.vx/v0 + 1.0)/self.vx

analyticalSolutionVelocity = dict((ci,DarcyVelocity(v=darcyVelocity_0*((-1.)**ci),vx=darcyVelocity_x*((-1.)**ci))) for ci in range(nc))
analyticalSolutionParticleVelocity = dict((ci,ParticleVelocity(v=darcyVelocity_0*((-1.)**ci),vx=darcyVelocity_x*((-1.)**ci))) for ci in range(nc))
analyticalSolutionParticleTrajectory=dict((ci,ParticleVelocityTrajectory(v=darcyVelocity_0*((-1.)**ci),vx=darcyVelocity_x*((-1.)**ci))) for ci in range(nc))
#analytical velocity function
darcyVelocityEval = dict((ci,lambda x,t,cii=ci: analyticalSolutionVelocity[cii].uOfXT(x,t)) for ci in range(nc))
#analytical solution
class GaussianSolution:
    def __init__(self,sigma=1.0/16.,xc=0.2,v=1.0,D=0.0):
        self.sigma= sigma
        self.xc   = xc
        self.v    = v
        self.D    = D
    def uOfXT(self,x,t):
        B = x[0] - self.v*t - self.xc
        if t <= 1.0e-12:
            return exp((-B**2)/(2.0*self.sigma**2))
        C  = sqrt(2.*self.sigma**2 + 4.0*self.D*t)
        A  = sqrt(2.0)*self.sigma/C
        return A*exp(-B**2/C**2)
#'reaction' term to enforce analytical solution (note it's the minus sign for the source)
class GaussianRHS:
    def __init__(self,sigma={0:1.0/16.},xc={0:0.2},v={0:1.0},D={0:0.0}):
        self.sigma= sigma
        self.xc   = xc
        self.v    = v
        self.D    = D
    def uOfXT_x(self,x,t,ci=0):
        B = x - self.v[ci]*t - self.xc[ci]
        C  = sqrt(2.*self.sigma[ci]**2 + 4.0*self.D[ci]*t)
        A  = sqrt(2.0)*self.sigma[ci]/C
        emB2dC2= numpy.exp(-B**2/C**2)
        C2 = C**2; C3 = C2*C; C4 = C3*C
        dudt = -2.0*sqrt(2.0)*self.sigma[ci]/C3*self.D[ci]*emB2dC2 + A*(2.0*B/C2*self.v[ci] + 4.0*self.D[ci]*B*B/C4)*emB2dC2
        advterm = -self.v[ci]*A*2.0*B/C2*emB2dC2
        dispterm= -self.D[ci]*A*emB2dC2*(-2.0/C2 + 4.0*B*B/C4)
        return dudt + advterm + dispterm
    def uOfXT(self,x,t,ci=0):
        return self.uOfXT_x(x[0],t,ci)
    def evaluate(self,t,c):
        #mwf debug
        #import pdb
        #pdb.set_trace()
        npoints = len(c['x'].flat)
        xind    = range(0,npoints,3)
        for ci in range(len(self.sigma)):
            c[('r',ci)].flat[:] = self.uOfXT_x(c['x'].flat[xind],t,ci)
#initial conditions 
sigma=dict((ci,0.0316) for ci in range(nc)); xc=dict((ci,0.3) for ci in range(nc))
for ci in range(nc):
    if darcyVelocity_0*((-1.)**ci) < 0.0:
        xc[ci] = 0.7
D = dict((ci,abs(velocity_0)*alpha_L + d) for ci in range(nc))
velocity_0_ci = dict((ci,velocity_0*((-1.)**ci)) for ci in range(nc))

coefficients = SubsurfaceTransportCoefficients.GroundwaterTransportCoefficientsELLAM(nc=nc,
                                                                                     omega=omega,
                                                                                     alpha_L=alpha_L,
                                                                                     alpha_T=alpha_T,
                                                                                     d = d,
                                                                                     nd = nd,
                                                                                     velocityFunctions = darcyVelocityEval,
                                                                                     velocitySpaceFlag = velocitySpaceFlag,
                                                                                     reactionTerms = GaussianRHS(sigma=sigma,xc=xc,v=velocity_0_ci,D=D),
                                                                                     flowIsTransient=velocityIsTransient)


smoothInflow = False; 
smoothDt=0.1; inflowVal=0.0
def inflowFlux(t,ci=0):
    if (smoothInflow and t <= smoothDt):
        return inflowVal*((-1.0)**(ci+1))*(t-0.0)/smoothDt
    return inflowVal*((-1.0)**(ci+1))

    

def getDBC(x,tag):
    pass
dirichletConditions = dict((ci,getDBC) for ci in range(nc))


#correct if constant velocity
analyticalSolution = dict((ci,GaussianSolution(sigma=sigma[ci],xc=xc[ci],v=velocity_0_ci[ci],D=D[ci])) for ci in range(nc))

initialConditions  = dict((ci,GaussianSolution(sigma=sigma[ci],xc=xc[ci],v=velocity_0_ci[ci],D=D[ci])) for ci in range(nc))

fluxBoundaryConditions = dict((ci,'outFlow') for ci in range(nc))#'noFlow'}

def getAFBC_ci(x,tag,ci=0):
    if x[0] == 0.0 and analyticalSolutionVelocity[ci].uOfXT(x,0.)[0] > 0.:
        return lambda x,t: inflowFlux(t,ci)
    if x[0] == 1.0 and analyticalSolutionVelocity[ci].uOfXT(x,0.)[0] < 0.:
        return lambda x,t: inflowFlux(t,ci)
advectiveFluxBoundaryConditions =  dict((cj,lambda x,tag,cii=cj : getAFBC_ci(x,tag,ci=cii)) for cj in range(nc))


def getDFBC_ci(x,tag,ci=0):
   if x[0] == 1.0 and analyticalSolutionVelocity[ci].uOfXT(x,0.)[0] > 0.:
       return lambda x,t: 0.0
   if x[0] == 0.0 and analyticalSolutionVelocity[ci].uOfXT(x,0.)[0] < 0.:
       return lambda x,t: 0.0
diffusiveFluxBoundaryConditions = {}
for ci in range(nc):
    diffusiveFluxBoundaryConditions[ci] = dict((ck,lambda x,tag,cii=ck : getDFBC_ci(x,tag,ci=cii)) for ck in range(nc))

T = 2.0*tForVelocityReversal#0.5#0.75
#tnList = [0,0.5]
#tnList = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
