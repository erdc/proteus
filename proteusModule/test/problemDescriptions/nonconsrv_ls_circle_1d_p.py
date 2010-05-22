from pyadh import *
from pyadh.default_p import *
from math import *
"""
Non-conservative advection of a circular (v-shaped) level set in 1D.
"""
##\page Tests Test Problems 
# \ref nonconsrv_ls_circle_1d_p.py "Non-conservative linear advection of v-shaped level set"
#

##\ingroup test
#\file nonconsrv_ls_circle_1d_p.py
#\brief Non-conservative advection of a circular (v-shaped) level set in 1D.
#
#The equation and coefficients are implemented in the ConstantVelocityLevelSet class.
#
nd = 1



class ConstantVelocityCircle:
    def __init__(self,radius=0.1,b=[1.],startX=0.25):
        self.radius = radius
        self.b      = b
        self.startX = startX
    def uOfXT(self,x,t):
        centerX = self.b[0]*t + self.startX
        return sqrt((x[0]-centerX)**2) - self.radius

b0 = Numeric.array([-1.],Numeric.Float)
x0 = 0.5
r0 = 1./8.
analyticalSolution = {0:ConstantVelocityCircle(r0,b0,x0)}

class ConstantVelocityLevelSet(TransportCoefficients.TC_base):
    """
    Coefficients for the non-conservative form of the level set equation with a constant velocity.
    """
    ##
    #\f{eqnarray*}
    #\phi_t + \mathbf{v} \cdot \grad \phi = 0
    #\f}
    def __init__(self,b=[1.0],lsModelId=0):
        self.b=b
        mass={0:{0:'linear'}}
#mwf need for grad(u)
        advection={0:{0:'linear'}}
        diffusion={}
        potential={}
        reaction={}
        hamiltonian={0:{0:'linear'}}
        TransportCoefficients.TC_base.__init__(self,
                                             1,
                                             mass,
                                             advection,
                                             diffusion,
                                             potential,
                                             reaction,
                                             hamiltonian)
        self.lsModelId = lsModelId
        self.lsModel= None  #the level set model itself
        self.tLast = -12345.0  #the last time step seen
    def attachModels(self,modelList):
        if len(modelList) > 1:
            self.lsModel = modelList[self.lsModelId]
    
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

        #mwf debug
        #print """L_S_CS_1D_circle evaluate t= %g b=%s """ % (t,self.b)
        for i in range(len(c[('u',0)].flat)):
            #print """c[('grad(u)',0)].flat[%d]= %s \n""" % (i,c[('grad(u)',0)].flat[i])
            c[('H',0)].flat[i]=self.b[0]*c[('grad(u)',0)].flat[i]
            c[('dH',0,0)].flat[nd*i:nd*(i+1)]=self.b[:]
         #end for
    #end def
                                          

coefficients = ConstantVelocityLevelSet(b0)

#mwf hack for level set redistancing
redistance=True
def finalizeStepRD(c):
    #print """ In finalizeStepRD c.lsModel= %s""" % c.lsModel
    #update level set coefficients recursive

    if c.lsModel != None and redistance:
        #mwf debug
        #print """calling coefficients calculate ..."""
        c.lsModel.calculateCoefficients()
        c.lsModel.updateTimeHistory()
    #if
#def
if redistance == True:
    finalizeStep = finalizeStepRD

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

T = 5.e-1
