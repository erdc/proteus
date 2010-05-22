from pyadh import *
from pyadh.default_p import *

nd = 2

U0     = 1.0  #m/s
height = 1.0 #m
length2height = 10.#10.0
length = length2height*height
#
Re     = 100.0
rho0   = 1.0
nu     = height*U0/Re
gms2   = 9.8 #m/s^2
Fr     = U0/sqrt(height*gms2)
mu     = nu*rho0 #should be 0.001
g      = [0.0,-(U0**2)/(height*Fr**2)]
#g      = [0.0,0.0]
print "NS porous Obstacle U0= %s Re= %s Fr= %s g= %s" % (U0,Re,Fr,g) 
meanGrainsize = 0.05 #0.05 #m
obstaclePorosity = 0.8#0.45#0.8#

inflow = 1.0 #(U0/U0)
L = (length2height,height/height,1.0)

obstacleXM = 0.45*length2height
obstacleXP = 0.55*length2height
obstacleYM = 0.0
obstacleYP = 0.4

def onSquareObstacle(x):
    if obstacleXM <= x[0] and x[0] <= obstacleXP:
        if obstacleYM <= x[1] and x[1] <= obstacleYP:
            return True
    return False

#trying LWI dike
#obstacleXM = 0.55*length2height
#obstacleXP = 0.6*length2height
slopeM = 1.0/6.0 #front slope
slopeP = 1.0/3.0 #back slope
obstacleDXM = (obstacleYP-obstacleYM)/slopeM
obstacleDXP = (obstacleYP-obstacleYM)/slopeP

def onLWIdike(x):
    if obstacleXM <= x[0] and x[0] <= obstacleXP:
        if obstacleYM <= x[1] and x[1] <= obstacleYP:
            return True
    if obstacleXM-obstacleDXM <= x[0] and x[0] <= obstacleXM:
        if obstacleYM <= x[1] and x[1] <= obstacleYM+slopeM*(x[0]-(obstacleXM-obstacleDXM)):
            return True
    if obstacleXP <= x[0] and x[0] <= obstacleXP+obstacleDXP:
        if obstacleYM <= x[1] and x[1] <= obstacleYP+slopeP*(obstacleXP-x[0]):
            return True
    return False

def setObstaclePorosity(x,porosity,meanGrain=None):
    porosity.flat[:]  = 1.0 #fluid domain
    if meanGrain != None:
        meanGrain.flat[:] = meanGrainsize
    #import pdb
    #pdb.set_trace()
    nPoints = len(x.flat)/3 #points always 3d
    for k in range(nPoints):
        if onSquareObstacle(x.flat[3*k:3*(k+1)]):
        #if onLWIdike(x.flat[3*k:3*(k+1)]):
            porosity.flat[k]=obstaclePorosity
            #
        #
    #
    #mwf hack
    #porosity.flat[:] = obstaclePorosity
#

coefficients = VolumeAveragedNavierStokesFullDevStress(rho=rho0,mu=mu,g=g,nd=nd,
                                                       meanGrainSize=meanGrainsize,
                                                       setParamsFunc=setObstaclePorosity,
                                                       stokesOnly=False)

#coefficients = NavierStokes(rho=rho0,nu=nu,g=g,nd=nd)
#coefficients = Stokes(rho=rho0,nu=nu,g=g,nd=nd)

T = 5.1e-2#5.1e-1#5.1# 5.1e-1
EPS=1.0e-6
def getDBC_p(x,tag):
    #works for all three cases
    #if x[0] == 0.0 and x[1] == 0.0:
    #    return lambda x,t: 0.0
    if x[0] >= L[0]-EPS:
        return lambda x,t: 0.0 + g[1]*rho0*(x[1]-0.0)
    #pass
    if x[1] >= L[1]-EPS:
        return lambda x,t: 0.0 + g[1]*rho0*(x[1]-0.0)
def getDBC_u(x,tag):
    if x[1] <= EPS:
        return lambda x,t: 0.0
    if x[0] <= EPS:
        return lambda x,t: inflow*(x[1]-0.0)/L[1]
    if x[1] >= L[1]-EPS:
        return lambda x,t: inflow
    
def getDBC_v(x,tag):
    if x[0] <= EPS:
        return lambda x,t: 0.0
    if x[1] <= EPS:
        return lambda x,t: 0.0
    #if x[1] == L[1]:
    #    return lambda x,t: 0.0
def getAFBC_p(x,tag):
    if x[1] <= EPS:# or x[1] == L[1]:
        return lambda x,t: 0.0
    if x[0] <= EPS:
        return lambda x,t: -inflow*(x[1]-0.0)/L[1]

def getDFBC_u(x,tag):
    if x[0] >= L[0]-EPS:
        return lambda x,t: 0.0
def getDFBC_v(x,tag):
    if x[0] >= L[0]-EPS:
        return lambda x,t: 0.0
#
def getNoBC(x,tag):
    pass


dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v}
    
#has nonzero v component
#fluxBoundaryConditions = {0:'outFlow',
#                          1:'mixedFlow',
#                          2:'outFlow'}
#ok too
#fluxBoundaryConditions = {0:'outFlow',
#                          1:'setFlow',
#                          2:'outFlow'}
#ok too
fluxBoundaryConditions = {0:'outFlow',
                          1:'outFlow',
                          2:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getNoBC,
                                    2:getNoBC}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{},#{1:getDFBC_u},
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
