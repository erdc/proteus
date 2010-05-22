from pyadh import *
from pyadh.default_p import *

import pome_square_pack_gen_poly

nd = 2

polyfile="pome_square"

#array of circles (cylinders)
ncircles=[4,4]

inflow = 0.1

periodic = False#True

T = 1000.0/inflow

lengthScale = 0.5/max(ncircles)

Re = 0.01 #inflow*L/nu

nu = abs(inflow)*lengthScale/Re

rho=1.0

coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                             sigma=0.0,
                                             rho_0=rho,nu_0=nu,
                                             rho_1=rho,nu_1=nu,
                                             g=[0.0,0.0,0.0],
                                             nd=nd,
                                             LS_model=None,
                                             KN_model=None,
                                             epsFact_density=None,
                                             stokes=False)

points_on_grain=21
points_on_boundary=41
dx = 1.0/float(points_on_boundary+1)
boundaryFlags = pome_square_pack_gen_poly.genPoly(polyfile,ncircles[0],ncircles[1],points_on_grain,points_on_boundary)
grainBoundaries = range(5,len(boundaryFlags)+1)

def getDBC_pressure_tube(x,flag):
    if flag == boundaryFlags['right']:
        return lambda x,t: (1.0-x[1])*coefficients.rho*coefficients.g[1]

def getDBC_u_tube(x,flag):
    if flag == boundaryFlags['left']:
        return lambda x,t: inflow
    if flag in grainBoundaries:
        return lambda x,t: 0.0
    if periodic:
        if flag in [boundaryFlags['back'],boundaryFlags['front']]:
            return lambda x,t: 0.0

def getDBC_v_tube(x,flag):
    if flag == boundaryFlags['left']:
        return lambda x,t: 0.0
    if flag in grainBoundaries:
        return lambda x,t: 0.0
    if periodic:
        if flag in [boundaryFlags['back'],boundaryFlags['front']]:
            return lambda x,t: 0.0

dirichletConditions = {0:getDBC_pressure_tube,
                       1:getDBC_u_tube,
                       2:getDBC_v_tube}


if periodic:
    def getPDBC(x,flag):
        if (flag == boundaryFlags['front'] or
            flag == boundaryFlags['back']):
            #make sure it's not a node on left or right boundary
            if x[0] > 0.0 and x[0] < 1.0:
                return numpy.array([x[0],0.0,x[2]])
    periodicDirichletConditions = {0:getPDBC,
                                   1:getPDBC,
                                   2:getPDBC}

def getAFBC_p_tube(x,flag):
    if flag in [boundaryFlags['back'],boundaryFlags['front']]:
        return lambda x,t: 0.0
    if flag == boundaryFlags['left']:
        return lambda x,t: -inflow
    if flag in grainBoundaries:
        return lambda x,t: 0.0

def getAFBC_vel_tube(x,flag):
    if periodic:
        if flag in [boundaryFlags['back'],boundaryFlags['front']]:
            return lambda x,t: 0.0

def getDFBC_vel_tube(x,flag):
    if periodic:
        if flag in [boundaryFlags['back'],boundaryFlags['front']]:
            return lambda x,t: 0.0

fluxBoundaryConditions = {0:'mixedFlow',
                          1:'mixedFlow',
                          2:'mixedFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_tube,
                                    1:getAFBC_vel_tube,
                                    2:getAFBC_vel_tube}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_vel_tube},
                                   2:{2:getDFBC_vel_tube}}

class Steady_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return (L[1]-x[1])*coefficients.rho*coefficients.g[1]

class Steady_u:
    def __init__(self,val=0.0):
        self.val=val
    def uOfXT(self,x,t):
        return self.val

class Steady_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v()}
