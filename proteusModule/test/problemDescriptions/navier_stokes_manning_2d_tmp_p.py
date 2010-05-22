from pyadh import *
from pyadh.default_p import *

import math


nd = 2
Size=[10.0,1.0]
angle=math.pi/10000.0#100.0
numarcs=5
amplitude=0.1




ncircles=[5,5]

Re=0.01


#inflow = 0.0 #0.1
#pInflow=100000.0
#pOutflow=0.0


rho=998.2
nu=1.004e-6

mu=rho*nu

T = 50.0

#lengthScale = 0.5/max(ncircles)

#Re = 0.01 #inflow*L/nu

#nu = 1.004e-6#abs(inflow)*lengthScale/Re

#rho=998.2
g=[9.8*math.cos(math.pi/2.0-angle),-9.8*math.sin(math.pi/2.0-angle),0.0]
#g=[0.0,-9.8,0.0]
coefficients = TwophaseNavierStokes_ST_LS_SO(epsFact=0.0,
                                             sigma=0.0,
                                             rho_0=rho,nu_0=nu,
                                             rho_1=rho,nu_1=nu,
                                             g=g,
                                             nd=nd,
                                             LS_model=None,
                                             KN_model=None,
                                             epsFact_density=None,
                                             stokes=False)
                                                   


points_on_boundary=7


## from manning2dDomain import*
## domain = manning2D( Lx=Size[0],
##                     Ly = Size[1],
##                     numarcs=numarcs,
##                     amplitude=amplitude,
##                     points_on_boundary=points_on_boundary)



from box2dDomain import *
domain = box2D(Lx=Size[0],Ly=Size[1])

domain.writePoly("box2D")
    
boundaryFlags=domain.boundaryFlags
grainBoundaries = range(5,len(boundaryFlags)+1)

dx = Size[0]/float(2*numarcs)
#boundaryFlags = pome_square_pack_gen_poly.genPoly(polyfile,ncircles[0],ncircles[1],points_on_grain,points_on_boundary)
#grainBoundaries = range(5,len(boundaryFlags)+1)


def getDBC_pressure_tube(x,flag):
    if flag == boundaryFlags['left']:
        #print 'left ' , x
        return lambda x,t: -rho*(Size[1]-x[1])*coefficients.g[1]
    if flag == boundaryFlags['right']:
        #print 'right ' , x
        return lambda x,t: -rho*(Size[1]-x[1])*coefficients.g[1]
    

def getDBC_u_tube(x,flag):
    #if flag in [boundaryFlags['back'],boundaryFlags['front']]:
    if flag == boundaryFlags['front']:
        return lambda x,t: 0.0
    #elif flag == boundaryFlags['left']:
    #    return lambda x,t: 0.0
##     if flag == boundaryFlags['left']:
##         return lambda x,t: x[1]/mu *math.tan(angle)*(Size[1]-x[1]/2.0)

def getDBC_v_tube(x,flag):
    #if flag in [boundaryFlags['front'],boundaryFlags['back']]:
    if flag == boundaryFlags['front']:
        return lambda x,t: 0.0
   
dirichletConditions = {0:getDBC_pressure_tube,
                       1:getDBC_u_tube,
                       2:getDBC_v_tube}



def getAFBC_p_tube(x,flag):
    if flag in [boundaryFlags['back'],boundaryFlags['front']]:
        return lambda x,t: 0.0

def getAFBC_vel_tube(x,flag):
    #pass
    #if flag in [boundaryFlags['back'],boundaryFlags['left'],boundaryFlags['right']]:
    #if flag in [boundaryFlags['back'],boundaryFlags['front']]:
     #   return lambda x,t: 0.0
     if flag == boundaryFlags['back']:
         return lambda x,t: 0.0
def getDFBC_vel_tube(x,flag):
    #pass
    if flag in [boundaryFlags['back'],boundaryFlags['left'],boundaryFlags['right']]:
    #if flag == boundaryFlags['back']:
        return lambda x,t: 0.0

fluxBoundaryConditions = {0:'setFlow',
                          1:'setFlow',
                          2:'setFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_p_tube,
                                    1:getAFBC_vel_tube,
                                    2:getAFBC_vel_tube}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_vel_tube,2:getDFBC_vel_tube},
                                   2:{1:getDFBC_vel_tube,2:getDFBC_vel_tube}}

class Steady_p:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return -rho*(Size[1]-x[1])*coefficients.g[1]
   
        

class Steady_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

class Steady_v:
    def __init__(self,val=0.0):
        pass
    def uOfXT(self,x,t):
        return 0.0


initialConditions = {0:Steady_p(),
                     1:Steady_u(),
                     2:Steady_v()}

class Couette_u:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return x[1]/mu *math.tan(angle)*(Size[1]-x[1]/2.0)
class Couette_v:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

## initialConditions = {0:Steady_p(),
##                      1:Couette_u(),
##                      2:Couette_v()}

## analyticalSolution= {0:Steady_p(),
##                      1:Couette_u(),
##                      2:Couette_v()}
                     
        
