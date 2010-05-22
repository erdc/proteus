from pyadh import *
from pyadh.default_p import *

nd = 2

height = 1.0 #m
length2height = 10.#10.0
length = length2height*height
U0     = 1.0 #m/s
gms2   = 9.8
Fr     = U0/sqrt(height*gms2)#[1 m/s]/[1 m]/[9.8 m/s^2]?
rho0   = 1.0
Re     = 1.0
nu     = U0*height/Re
mu     = nu*rho0
g      = [0.0,-height/(gms2*Fr**2)]
#g      = [0.0,0.0]
print "S porous Obstacle Re= %s Fr= %s g= %s" % (Re,Fr,g) 
meanGrainsize = 0.05 #m
obstaclePorosity = 0.45#0.5
#alphaP        = 0.0#1.e3
#betaP         = 0.0#1.1
#gammaP        = 0.34
#KeuleganCarpenter = 1000.0#Keulegan Carpenter Number u_c(U0?)*T/(porosity*meanGrainsize)

inflow = 1.0 #(U0/U0)
L = (length2height,height/height,1.0)

obstacleXM = 0.45*length2height
obstacleXP = 0.55*length2height
obstacleYM = 0.0
obstacleYP = 0.4*height

def setObstaclePorosity(x,porosity):
    porosity.flat[:]  = 1.0 #fluid domain
    #import pdb
    #pdb.set_trace()
    nPoints = len(x.flat)/3 #points always 3d
    for k in range(nPoints):
        if obstacleXM <= x.flat[k*3+0] and x.flat[k*3+0] <= obstacleXP:
            if obstacleYM <= x.flat[k*3+1] and x.flat[k*3+1] <= obstacleYP:
                porosity.flat[k]=obstaclePorosity
            #
        #
    #
    #mwf hack
    #porosity.flat[:] = obstaclePorosity
#

#coefficients = VolumeAveragedNavierStokes(rho=rho0,nu=nu,g=g,nd=nd,
#                                          meanGrainSize=meanGrainsize,
#                                          alphaP=alphaP,
#                                          betaP=betaP,
#                                          gammaP=gammaP,
#                                          KeuleganCarpenter=KeuleganCarpenter,
#                                          setParamsFunc=setObstaclePorosity,
#                                          stokesOnly=True)
coefficients = VolumeAveragedNavierStokesFullDevStress(rho=rho0,mu=mu,g=g,nd=nd,
                                                       meanGrainSize=meanGrainsize,
                                                       setParamsFunc=setObstaclePorosity,
                                                       stokesOnly=True)

#coefficients = Stokes(rho=rho0,nu=nu,g=g,nd=nd)

T = 1.0

def getDBC_p(x):
    #works for all three cases
    #if x[0] == 0.0 and x[1] == 0.0:
    #    return lambda x,t: 0.0
    if x[0] == L[0]:
        return lambda x,t: 0.0 + g[1]*rho0*(x[1]-0.0)
    #pass
    if abs(x[1]-L[1]) < 1.0e-6:
        return lambda x,t: 0.0 + g[1]*rho0*(x[1]-0.0)
def getDBC_u(x):
    if x[1] == 0.0:
        return lambda x,t: 0.0
    if x[0] == 0.0:
        return lambda x,t: inflow*(x[1]-0.0)/L[1]
    if x[1] == L[1]:
        return lambda x,t: inflow
    
def getDBC_v(x):
    if x[0] == 0.0:
        return lambda x,t: 0.0
    if fabs(x[1]-0.0) < 1.0e-4:
        return lambda x,t: 0.0
    #if x[1] == L[1]:
    #    return lambda x,t: 0.0
    #force to be zero on outflow
    #if x[0] == L[0]:
    #    return lambda x,t: 0.0
def getAFBC_p(x):
    if x[1] == 0.0:# or x[1] == L[1]:
        return lambda x,t: 0.0
    if x[0] == 0.0:
        return lambda x,t: -inflow*(x[1]-0.0)/L[1]

def getDFBC_u(x):
    if x[0] == L[0]:
        return lambda x,t: 0.0
def getDFBC_v(x):
    if x[0] == L[0]:
        return lambda x,t: 0.0
#
def getNoBC(x):
    pass


dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v}
    
#has nonzero v component
#fluxBoundaryConditions = {0:'outFlow',
#                          1:'mixedFlow',
#                          2:'outFlow'}
#ok too
fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getNoBC,
                                    2:getNoBC}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{}}#or 2:getDFBC_v same behavior

class Hydrostatic_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0 + g[1]*rho0*(x[1]-0.0)

class NoFlow:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Couette_u:
    def __init__(self,U0=inflow):
        self.U0=inflow
    def uOfXT(self,x,t):
        return self.U0*(x[1]-0.0)/L[1]

initialConditions = {0:Hydrostatic_p(),
                     1:Couette_u(),
                     2:NoFlow()}
