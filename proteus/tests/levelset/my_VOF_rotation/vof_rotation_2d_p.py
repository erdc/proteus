from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from rotation2D import *
from proteus.mprans import VOF
name=soname+"_vof"

"""
The non-conservative level set description of a bubble in a two-phase flow
"""

LevelModelType = VOF.LevelModel

#My Own Coefficients
class MyCoefficients(VOF.Coefficients):
    def attachModels(self,modelList):
        self.model = modelList[self.modelIndex]
	self.u_old_dof = numpy.copy(self.model.u[0].dof)
        self.q_v = numpy.zeros((self.model.mesh.nElements_global,self.model.nQuadraturePoints_element,self.model.nSpace_global),'d')    
        self.ebqe_v = numpy.zeros((self.model.mesh.nExteriorElementBoundaries_global,self.model.nElementBoundaryQuadraturePoints_elementBoundary),'d')
        self.model.q[('velocity',0)]=self.q_v
        self.model.ebqe[('velocity',0)]=self.ebqe_v
        self.ebqe_phi = numpy.zeros(self.model.ebqe[('u',0)].shape,'d') #NOTE: this is not needed (is for LS)
    def preStep(self,t,firstStep=False):
        pi = math.pi
        import numpy as np
        x = self.model.q['x'][...,0]
        y = self.model.q['x'][...,1]
        #ROTATION
        self.q_v[...,0]  = -2.0*pi*y
        self.q_v[...,1]  =  2.0*pi*x
        #PERIODIC VORTEX
        #T=8
        #self.q_v[...,0] = -2*np.sin(pi*y)*np.cos(pi*y)*np.sin(pi*x)**2*np.cos(pi*t/T)
        #self.q_v[...,1] = 2*np.sin(pi*x)*np.cos(pi*x)*np.sin(pi*y)**2*np.cos(pi*t/T)        
        #TRANSLATION
        #self.q_v[...,0]  = 0.0
        #self.q_v[...,1]  = -1.0
        copyInstructions = {}
        return copyInstructions
    def postStep(self,t,firstStep=False):
       	self.u_old_dof = numpy.copy(self.model.u[0].dof)
        copyInstructions = {}
        return copyInstructions
    def evaluate(self,t,c):
        pass

#coefficients = VOF.Coefficients(RD_model=None,ME_model=0,checkMass=checkMass,
#                                    epsFact=epsFact_vof,useMetrics=useMetrics)
coefficients = MyCoefficients(epsFact=epsFactHeaviside,checkMass=checkMass,useMetrics=useMetrics,ME_model=0,
                              cE=cE,cMax=cMax,cK=cK,ENTROPY_VISCOSITY=ENTROPY_VISCOSITY,SUPG=SUPG)
#coefficients = MyCoefficients(epsFact=epsFactHeaviside,checkMass=checkMass,useMetrics=useMetrics,ME_model=0)

def Heaviside(phi):
    if phi > 0:
        return 1.0
    elif phi < 0:
        return 0.0
    else:
        return 0.5

class init_cond:
    def __init__(self,center=[0.5,0.75,0.5],radius=0.15):
        self.radius  = radius
        self.center  = center
    def uOfX(self,X):
        dx = X[0]-self.center[0]; dy = X[1]-self.center[1];
        dBubble = self.radius - sqrt(dx**2 + dy**2)

        #Zalesak disk
        #dBubble = 1.0*(sqrt(dx**2 + dy**2) <= self.radius) - 1.0*(sqrt(dx**2 + dy**2) > self.radius)
        #xSlit1 = X[0] < self.center[0]+0.025
        #xSlit2 = X[0] > self.center[0]-0.025
        #xSlit = xSlit1 and xSlit2
        #ySlit = X[1] < 0.75+0.1125
        #slit = xSlit*ySlit
        #if (slit==1):
        #    dBubble = -1
        return smoothedHeaviside(epsFactHeaviside*he,dBubble)#Heaviside(dBubble)
    #end
    def uOfXT(self,X,t):
        return self.uOfX(X)
    #end
#end init_cond

analyticalSolutions = None

def getDBC(x,flag):
    pass

dirichletConditions = {0:getDBC}

initialConditions  = {0:init_cond(center=[0.0,0.5],radius=0.25)}
#initialConditions  = {0:init_cond(center=[0.5,0.75],radius=0.15)}

fluxBoundaryConditions = {0:'outFlow'}

#cek made no flux since v.n = 0 for this v
def getAFBC(x,flag):
   return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:getAFBC}

diffusiveFluxBoundaryConditions = {0:{}}
