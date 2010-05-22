from pyadh import *
from pyadh.default_p import *

nd = 2

U0     = 1.0  #m/s
height = 1.2  #m
length = 51.0 #m
length2height = length/height#10.0
#
Re     = 10.0#100.0
rho0   = 1.0
nu     = height*U0/Re
gms2   = 9.8 #m/s^2
Fr     = U0/sqrt(height*gms2)
mu     = nu*rho0 #should be 0.001
g      = [0.0,-(U0**2)/(height*Fr**2)]
#g      = [0.0,0.0]
print "NS porous Obstacle U0= %s Re= %s Fr= %s g= %s" % (U0,Re,Fr,g) 
meanGrainsize = 0.05 #0.05 #m
obstaclePorosity = 0.45#0.45#0.8#

inflow = 1.0 #(U0/U0)
L = (length,height,1.0)


polyfile = 'cobras_wave_flume_1'

porosityTypes = {0 : 1.0, #fluid
                 1 : 0.5, #tetrapod
                 2 : 0.53,#filter_layer
                 3 : 0.49,#rubble_mound
                 4 : 1.0e-8}#caisson not meshed
meanGrainTypes= {0 : 0.05,#fluid shouldn't matter
                 1 : 0.05, #tetrapod
                 2 : 0.05,#filter_layer
                 3 : 0.05,#rubble_mound
                 4 : 0.05}#caisson not meshed
                 
coefficients = VolumeAveragedNavierStokesFullDevStress(rho=rho0,mu=mu,g=g,nd=nd,
                                                       meanGrainSize=meanGrainsize,
                                                       setParamsFunc=None,
                                                       stokesOnly=False,
                                                       meanGrainSizeTypes=meanGrainTypes,
                                                       porosityTypes=porosityTypes)

#coefficients = NavierStokes(rho=rho0,nu=nu,g=g,nd=nd)
#coefficients = Stokes(rho=rho0,nu=nu,g=g,nd=nd)

T = 5.1e-2#5.1# 5.1e-1
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
