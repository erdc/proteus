from pyadh import *
from pyadh.default_p import *
"""
3D,Linear advection of a Gaussian with ELLAM
"""

## \page Tests Test Problems 
# \ref la_gauss_ellam_3d_p.py "Linear advection-diffusion of a Gaussian with ELLAM"
#

##\ingroup test
#\file la_gauss_ellam_3d_p.py
#@{
#
# \brief Uncoupled, multicomponent linear advection-diffusion for spherical initial conditions.
#   If rotating velocity, clockwise for ci % 2 = 0, counterclockwise for ci % 2 = 1
#   If not rotating velocity, scale basic velocity by (-1)**ci
# \todo finish la_gauss_3d_p.py doc
from pyadh import LADRellam
LevelModelType = LADRellam.OneLevelLADR
name = "la_gauss_ellam_3d"

nd = 3
nc = 1
#dispersivities
alpha_L = 1.0e-4
alpha_T = 1.0e-5
#molecular diffusion
d       = 0.0
#porosity
omega   = 1.0

#constant and possible derivative in space for velocities
darcyVelocity_0 = numpy.array([0.0,1.0,1.0])
darcyVelocity_x = 0.0
velocity_0      = darcyVelocity_0/omega
velocity_x      = darcyVelocity_x/omega
rotatingVelocity   = False
velocityIsTransient= False
tForVelocityReversal = 0.75#10000.0
L = (1.0,1.0,1.0)
modelColumn = False
if rotatingVelocity == False:
    modelColumn = False#True
if modelColumn:
    L = (0.1,0.1,1.0)
if rotatingVelocity:
    class DarcyVelocity:
        """
        rotating velocity in x,y constant in z
        """
        def __init__(self,xc=0.5,yc=0.5,zc=0.5,zVelocity=0.0,
                     transient=velocityIsTransient,tForReversal=tForVelocityReversal,
                     clock=1.0):
            self.xc=xc; self.yc=yc; self.zc=zc; self.zVelocity=0.0; self.clock=clock
            self.transient=transient; self.tForReversal=tForReversal
        def uOfXT(self,x,t):
            #mwf debug
            #print "la_gauss_ellam_3d_p velocity v= %s vx= %s x= %s returning %s \n" % (self.v,self.vx,x,self.v+self.vx*x[:nd])
            if self.transient:
                v = [2.0*pi*(x[1]-self.xc),2.0*pi*(self.yc-x[0]),self.zVelocity]
                return numpy.array(v)*(self.tForReversal-t)/(self.tForReversal-0.0)*self.clock
            return numpy.array([2.0*pi*(x[1]-self.xc),2.0*pi*(self.yc-x[0]),self.zVelocity])*self.clock
    class ParticleVelocity:
        def __init__(self,xc=0.5,yc=0.5,zc=0.5,zVelocity=0.0,transient=velocityIsTransient,tForReversal=tForVelocityReversal,
                     clock=1.0):
            self.xc=xc; self.yc=yc; self.zc=zc; self.zVelocity=zVelocity; self.clock = clock
            self.transient=transient; self.tForReversal=tForReversal
        def uOfXT(self,x,t):
            #mwf debug
            #print "la_gauss_ellam_3d_p velocity v= %s vx= %s x= %s returning %s \n" % (self.v,self.vx,x,self.v+self.vx*x[:nd])
            if self.transient:
                v = [2.0*pi*(x[1]-self.xc),2.0*pi*(self.yc-x[0]),self.zVelocity]
                return numpy.array(v)*(self.tForReversal-t)/(self.tForReversal-0.0)*self.clock
            return numpy.array([2.0*pi*(x[1]-self.xc),2.0*pi*(self.yc-x[0]),self.zVelocity])*self.clock
    class ParticleVelocityTrajectory:
        """
        not right, need to finish
        """
        def __init__(self,xc=0.5,yc=0.5,clock=1.0):
            self.xc=xc ; self.yc = yc
        def uOfXT(self,x,dt):
            x0 = x[0]-self.xc; y0 = x[1] - self.yc
            r0 = numpy.sqrt(x0**2 + y0**2)
            if abs(x0) < 1.e-12:
                if y0 >= 0.0:
                    th0 = 0.5*pi
                else:
                    th0 = 1.5*pi
            else:
                th0 = atan(y0/x0)
                if x0 < 0.0:
                    th0 = pi+th0
            th = th0 - 2.*pi*dt
            
            #mwf debug
            #print "la_gauss_ellam_3d_p velocityTrajectory x= %s dt= %s r0 = %s th0=%s returning x= %s \n" % (x,dt,r0,th0,
            #                                                                                                 numpy.array([r0*cos(th) + self.xc,r0*sin(th0) + self.yc]))
            return numpy.array([r0*cos(th) + self.xc,r0*sin(th) + self.yc,self.zVelocity])
        def dtOfdX(self,x,dx,n):
            """
            for now just take nearest possible intersection point, take advantage of rectangular domain
            """
            return 0.0
                 

else: # useConstantVelocity:
    class DarcyVelocity:
        def __init__(self,v=darcyVelocity_0,vx=darcyVelocity_x,transient=velocityIsTransient,tForReversal=tForVelocityReversal):
            self.v = v; self.vx = vx
            self.transient=transient; self.tForReversal=tForReversal
        def uOfXT(self,x,t):
            #mwf debug
            #print "la_gauss_ellam_3d_p velocity v= %s vx= %s x= %s returning %s \n" % (self.v,self.vx,x,self.v+self.vx*x[:nd])
            if self.transient:
                v = (self.v + self.vx*x[:nd])*(self.tForReversal-t)/(self.tForReversal-0.0)
                return v
            return self.v + self.vx*x[:nd]
    class ParticleVelocity:
        def __init__(self,v=velocity_0,vx=velocity_x,transient=velocityIsTransient,tForReversal=tForVelocityReversal):
            self.v = v; self.vx = vx
            self.transient=transient; self.tForReversal=tForReversal
        def uOfXT(self,x,t):
            #mwf debug
            #print "la_gauss_ellam_3d_p velocity v= %s vx= %s x= %s returning %s \n" % (self.v,self.vx,x,self.v+self.vx*x[:nd])
            if self.transient:
                v = (self.v + self.vx*x[:nd])*(self.tForReversal-t)/(self.tForReversal-0.0)
                return v
            return self.v + self.vx*x[:nd]
    class ParticleVelocityTrajectory:
        def __init__(self,v=velocity_0,vx=velocity_x):
            self.v = v; self.vx = vx
        def uOfXT(self,x,dt):
            if abs(self.vx) < 1.0e-10:
                return x[:nd] + dt*self.v
            else:
                v0 = self.v + self.vx*x[:nd]
                return x[:nd] + v0*(exp(self.vx*dt)-1.0)/self.vx
        def dtOfdX(self,x,dx,n):
            v0 = self.v + self.vx*x[:nd]
            dxm=numpy.dot(dx,n)
            v0m=numpy.dot(v0,n)
            if abs(v0m) < 1.0e-10:
                return 0.0
            if abs(self.vx) < 1.0e-10:
                return dxm/v0m
            if dxm*self.vx/v0m <= -1: return None
            return math.log(dxm*self.vx/v0m + 1.0)/self.vx

if rotatingVelocity == True:
    analyticalSolutionVelocity = dict((ci,DarcyVelocity(clock=(-1.0)**(ci))) for ci in range(nc))
    analyticalSolutionParticleVelocity = dict((ci,ParticleVelocity(clock=(-1.0)**(ci))) for ci in range(nc))
    analyticalSolutionParticleTrajectory=dict((ci,ParticleVelocityTrajectory(clock=(-1.0)**(ci))) for ci in range(nc))

else:
    analyticalSolutionVelocity = dict((ci,DarcyVelocity(v=darcyVelocity_0*((-1.)**ci),vx=darcyVelocity_x*((-1.)**ci))) for ci in range(nc))
    analyticalSolutionParticleVelocity = dict((ci,ParticleVelocity(v=darcyVelocity_0*((-1.)**ci),vx=darcyVelocity_x*((-1.)**ci))) for ci in range(nc))
    analyticalSolutionParticleTrajectory=dict((ci,ParticleVelocityTrajectory(v=darcyVelocity_0*((-1.)**ci),vx=darcyVelocity_x*((-1.)**ci))) for ci in range(nc))
#now define the Dirichlet boundary conditions
darcyVelocityEval = dict((ci,lambda x,t,cii=ci: analyticalSolutionVelocity[cii].uOfXT(x,t)) for ci in range(nc))

velocitySpaceFlag = 'c0p1'#'bdm1'#'bdm1'#'c0p1' #'rt0'

coefficients = SubsurfaceTransportCoefficients.GroundwaterTransportCoefficientsELLAM(nc=nc,
                                                                                     omega=omega,
                                                                                     alpha_L=alpha_L,
                                                                                     alpha_T=alpha_T,
                                                                                     d = d,
                                                                                     nd = nd,
                                                                                     velocityFunctions = darcyVelocityEval,
                                                                                     velocitySpaceFlag = velocitySpaceFlag,
                                                                                     flowIsTransient=velocityIsTransient)


smoothInflow = False; smoothDt=0.1; inflowVal=0.0
def inflowFlux(t,ci=0):
    if (smoothInflow and t <= smoothDt):
        return inflowVal*((-1.0)**(ci+1))*(t-0.0)/smoothDt
    return inflowVal*((-1.0)**(ci+1))

    

def getDBC(x,tag):
    pass
dirichletConditions = dict((ci,getDBC) for ci in range(nc))

class GaussIC:
    def __init__(self,sigma=1.0/16.,xc=[0.2,0.2,0.2],b=[0.0,0.0,1.0]):
        self.sigma= sigma
        self.xc   = numpy.array(xc)
        self.b    = numpy.array(b)
    def uOfXT(self,x,t):
        xct= self.xc + self.b*t
        d2 = numpy.sum((x[:nd]-xct)**2)
        return exp(-0.5*d2/self.sigma**2)

#initial conditions 
xc=dict((ci,[0.3,0.5,0.3]) for ci in range(nc))
if modelColumn:
    xci=dict((ci,[0.5*L[0],0.5*L[1],0.3*L[2]]) for ci in range(nc))
sigma=dict((ci,0.0316) for ci in range(nc))

for ci in range(nc):
    if not rotatingVelocity and darcyVelocity_0[2]*((-1.)**ci) < 0.0:
        xc[ci] = [0.3,0.5,0.7]

analyticalSolution = dict((ci,GaussIC(sigma=sigma[ci],xc=xc[ci])) for ci in range(nc))

initialConditions  = dict((ci,GaussIC(sigma=sigma[ci],xc=xc[ci])) for ci in range(nc))

fluxBoundaryConditions = dict((ci,'outFlow') for ci in range(nc))#'noFlow'}

def getAFBC_ci(x,tag,ci=0):
    if x[2] == 0.0 and analyticalSolutionVelocity[ci].uOfXT(x,0.)[2] > 0.:
        return lambda x,t: inflowFlux(t,ci)
    if x[2] == 1.0 and analyticalSolutionVelocity[ci].uOfXT(x,0.)[2] < 0.:
        return lambda x,t: inflowFlux(t,ci)
advectiveFluxBoundaryConditions =  dict((cj,lambda x,tag,cii=cj : getAFBC_ci(x,tag,ci=cii)) for cj in range(nc))

def getDFBC_ci(x,tag,ci=0):
   if x[2] == 1.0 and analyticalSolutionVelocity[ci].uOfXT(x,0.)[2] > 0.:
       return lambda x,t: 0.0
   if x[2] == 0.0 and analyticalSolutionVelocity[ci].uOfXT(x,0.)[2] < 0.:
       return lambda x,t: 0.0

diffusiveFluxBoundaryConditions = {}
for ci in range(nc):
    diffusiveFluxBoundaryConditions[ci] = dict((ck,lambda x,tag,cii=ck : getDFBC_ci(x,tag,ci=cii)) for ck in range(nc))


T = 0.75#0.75#2.0*0.75#0.75

