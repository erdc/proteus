from pyadh import *
from pyadh.default_p import *

dir=True

pInflow =1000.0

nd = 2


ellipse=True
b_over_a = 1.0/0.8
tilt_angle= 0.0

randomcylinder=False
num_spheres=25

#array of circles (cylinders)
ncircles=[5,5]

inflow = 0.0

periodic = True

T =1000000.0 #1000.0/inflow

lengthScale = 0.5/max(ncircles)

Re = 0.01 #inflow*L/nu

#nu = 1.0#abs(inflow)*lengthScale/Re

#rho=1.0
rho=1.0
nu=1.0 #1.0e-3
#rho_1=1.205,nu_1=1.500e-5,

if dir:
    g1=0.0
    g2=1.0
else:
    g1=1.0
    g2=0.0
    


coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                             sigma=0.0,
                                             rho_0=rho,nu_0=nu,
                                             rho_1=rho,nu_1=nu,
                                             g=[g1,g2,0.0],
                                             nd=nd,
                                             LS_model=None,
                                             KN_model=None,
                                             epsFact_density=None,
                                             stokes=True)

points_on_grain=30
points_on_boundary=50
dx = 1.0/float(points_on_boundary+1)
#boundaryFlags = pome_square_pack_gen_poly.genPoly(polyfile,ncircles[0],ncircles[1],points_on_grain,points_on_boundary)

if ellipse == True:
    from ellipse2dDomain import*
    domain = ellipse2D( nx= ncircles[0],
                     ny= ncircles[1],
                     b_over_a = b_over_a,
                     tilt_angle=tilt_angle,
                     points_on_grain=points_on_grain,
                     points_on_boundary= points_on_boundary)
    domain.writePoly("ellipse2D")
elif randomcylinder== True:
    from pome2dRandomDomain import*
    domain= pome2DRandom( num_spheres=num_spheres,
                          points_on_grain=points_on_grain,
                          points_on_boundary=points_on_boundary)
    domain.writePoly("pome2DRandom")
else:
    from pome2dDomain import*
    domain = pome2D( nx= ncircles[0],
                     ny= ncircles[1],
                     points_on_grain= points_on_grain,
                     points_on_boundary= points_on_boundary)

    domain.writePoly("pome2D")
    
boundaryFlags=domain.boundaryFlags
grainBoundaries = range(5,len(boundaryFlags)+1)
#pInflow =0.20* 1.0*0.1*998.2
pInflowback = 1.0*0.1*998.2
pOutflow = 0.0
def getDBC_pressure_tube(x,flag):
#    if flag == boundaryFlags['right']:
#        return lambda x,t: (1.0-x[1])*coefficients.rho*coefficients.g[1]
    if dir:
        if flag == boundaryFlags['front']:
            return lambda x,t: pInflow
        if flag == boundaryFlags['back']:
            return lambda x,t: pOutflow
    else:
        if flag == boundaryFlags['right']:
            return lambda x,t: pOutflow
        if flag == boundaryFlags['left']:
            return lambda x,t: pInflow
##    if flag== boundaryFlags['front']:
##        return lambda x,t: pInflowback
##    if flag== boundaryFlags['back']:
##        return lambda x,t: pOutflow
        
def getDBC_u_tube(x,flag):
#    if flag == boundaryFlags['left']:
#        return lambda x,t: inflow
    if flag in grainBoundaries:
        return lambda x,t: 0.0
    if periodic:
        if dir:
            if flag in [boundaryFlags['right'], boundaryFlags['left']]:
                return lambda x,t: 0.0
        else:
            if flag in [boundaryFlags['back'],boundaryFlags['front']]:
                return lambda x,t: 0.0
       

def getDBC_v_tube(x,flag):
#    if flag == boundaryFlags['left']:
#        return lambda x,t: 0.0
    if flag in grainBoundaries:
        return lambda x,t: 0.0
    if periodic:
        if dir:
            if flag in [boundaryFlags['right'], boundaryFlags['left']]:
                return lambda x,t: 0.0
        else:
            if flag in [boundaryFlags['back'],boundaryFlags['front']]:
                return lambda x,t: 0.0

dirichletConditions = {0:getDBC_pressure_tube,
                       1:getDBC_u_tube,
                       2:getDBC_v_tube}


if periodic:
    def getPDBC(x,flag):
        if dir:
            if (flag == boundaryFlags['left'] or
                flag == boundaryFlags['right']):
            #make sure it's not a node on left or right boundary
                if x[1] > 0.0 and x[1] < 1.0:
                    return numpy.array([0.0,x[1],x[2]])
        else:
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
 #   if flag == boundaryFlags['left']:
 #       return lambda x,t: -inflow
    if flag in grainBoundaries:
        return lambda x,t: 0.0

def getAFBC_vel_tube(x,flag):
    if periodic:
        if dir:
            if flag in [boundaryFlags['left'],boundaryFlags['right']]:
                return lambda x,t: 0.0
        else:
            if flag in [boundaryFlags['back'],boundaryFlags['front']]:
                return lambda x,t: 0.0

def getDFBC_vel_tube(x,flag):
    if periodic:
        if dir:
            if flag in [boundaryFlags['left'],boundaryFlags['right']]:
                return lambda x,t: 0.0
        else:
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
      # return (L[1]-x[1])*coefficients.rho*coefficients.g[1]
       return 0.0
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
