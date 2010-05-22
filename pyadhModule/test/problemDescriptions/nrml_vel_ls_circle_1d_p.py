from pyadh import *
from pyadh.default_p import *
from math import *
"""
Transport of a circular level set function by unit speed  normal motion.
"""
##\page Tests Test Problems 
# \ref nrml_vel_ls_circle_1d_p.py "Level set transport by unit speed normal motion (travel time equation)"
#

##\ingroup test
#\file nrml_vel_ls_circle_1d_p.py
#\brief Transport of a circular (v-shaped) level set function by unit speed  normal motion.
#
#The equation formulation and coefficients are implemented in the UnitNormalVelocityLevelSet class.
#
#\todo move travel time coefficients into TransportCoefficients.py

nd = 1



class UnitNormalVelocityCircle:
    def __init__(self,radius=0.1,startX=0.5):
        self.radius = radius
        self.startX = startX
    def uOfXT(self,x,t):
        centerX = self.startX
        return sqrt((x[0]-centerX)**2) - self.radius

x0 = 0.5
r0 = 1./8.
analyticalSolution = {0:UnitNormalVelocityCircle(r0,x0)}

class UnitNormalVelocityLevelSet(TransportCoefficients.TC_base):
    def __init__(self):
        mass={0:{0:'linear'}}
#mwf need for grad(u)
        advection={0:{0:'linear'}}
        diffusion={}
        potential={}
        reaction={}
        hamiltonian={0:{0:'nonlinear'}}
        TransportCoefficients.TC_base.__init__(self,
                                             1,
                                             mass,
                                             advection,
                                             diffusion,
                                             potential,
                                             reaction,
                                             hamiltonian)
    def evaluate(self,t,c):
        if ('m',0) in c.keys():
            c[('m',0)].flat[:] = c[('u',0)].flat[:]
        if ('dm',0,0) in c.keys():
            c[('dm',0,0)].flat[:] = 1.0
        if ('f',0) in c.keys():
            c[('f',0)].flat[:] = 0.0
        if ('df',0,0) in c.keys():
            c[('df',0,0)].flat[:] = 0.0
        if ('r',0) in c.keys():
            c[('r',0)].flat[:] = 0.0

        for i in range(len(c[('u',0)].flat)):
            agrad = sqrt(c[('grad(u)',0)].flat[i]**2)
            c[('H',0)].flat[i]=agrad
            c[('dH',0,0)].flat[nd*i:nd*(i+1)]=c[('grad(u)',0)].flat[i]/(agrad+1.0e-8)#sgrad*agrad
        #end for
    #end def
                                          

coefficients = UnitNormalVelocityLevelSet()

coefficients.variableNames=['u']

#now define the Dirichlet boundary conditions

def getDBC(x):
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

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 5.0e-1

