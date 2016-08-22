from proteus import *
from proteus.default_p import *
from proteus import Context
ct = Context.get()
from floating_bar import *
from proteus.mprans import MoveMesh

initialConditions = None

analyticalSolution = {}

nMediaTypes=1
smTypes      = numpy.zeros((nMediaTypes+1,2),'d')
smFlags      = numpy.zeros((nMediaTypes+1,),'i')

smTypes[0,0] = 1.0    ##E
smTypes[0,1] = 0.3    ##nu
smTypes[1,0] = 1.0    ##E
smTypes[1,1] = 0.3    ##nu

LevelModelType = MoveMesh.LevelModel
coefficients = MoveMesh.Coefficients(nd = ct.nd,
				     V_model=1,modelType_block=smFlags,
				     modelParams_block=smTypes,meIndex=0)
class FloatingObstacle(AuxiliaryVariables.AV_base):
    import numpy as np
    def __init__(self):
        self.object = None
    def attachModel(self,model,ar):
        self.model=model
        return self
    def attachAuxiliaryVariables(self,avDict):
        self.object = avDict['twp_navier_stokes_p'][0]
    def hx(self,x,t):
        if self.object == None:
            return 0.0
        else:
            hx = self.object.body.getRelPointPos(np.dot(self.object.last_rotation_inv,(x-self.object.last_position)))[0] - x[0]
#            hcx = self.object.h[0]
#            if fabs(hx-hcx)/(fabs(hcx)+1.0e-8) > 1.0e-8:
#                print "hx hcx",hx,hcx
            return hx
    def hy(self,x,t):
        if self.object == None:
            return 0.0
        else:
            hy = self.object.body.getRelPointPos(np.dot(self.object.last_rotation_inv,(x-self.object.last_position)))[1] - x[1]
#            hy = self.object.body.getPointVel(self.object.last_rotation_inv*(x-self.object.last_position))[1]*self.object.model.stepController.dt_model
#            hcy = self.object.h[1]
#            if fabs(hy-hcy)/(fabs(hcy)+1.0e-8) > 1.0e-8:
#                print "hy hcy",hy,hcy
            return hy
    def hz(self,x,t):
        if self.object == None:
            return 0.0
        else:
            hz = self.object.body.getRelPointPos(np.dot(self.object.last_rotation_inv,(x-self.object.last_position)))[2] - x[2]
            #print "z",x[2],np.dot(self.object.last_rotation_inv,(x-self.object.last_position))[2],hz
#            hz = self.object.body.getPointVel(self.object.last_rotation_inv*(x-self.object.last_position))[2]*self.object.model.stepController.dt_model
#            hcz = self.object.h[2]
#            if fabs(hz-hcz)/(fabs(hcz)+1.0e-8) > 1.0e-8:
#                print "hz hcz",hz,hcz
            return hz
    def calculate(self):
        pass

fo = FloatingObstacle()

def getDBC_hx(x,flag):
    if flag in [boundaryTags['left'],boundaryTags['right'], boundaryTags['top'], boundaryTags['bottom'], boundaryTags['front'], boundaryTags['back']]:
        return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        return lambda x,t: fo.hx(x,t)

def getDBC_hy(x,flag):
    #if flag in [boundaryTags['front'],boundaryTags['back']]:
    if flag in [boundaryTags['left'],boundaryTags['right'], boundaryTags['top'], boundaryTags['bottom'], boundaryTags['front'], boundaryTags['back']]:
        return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        return lambda x,t: fo.hy(x,t)

def getDBC_hz(x,flag):
    #if flag in [boundaryTags['top'],boundaryTags['bottom']]:
    if flag in [boundaryTags['left'],boundaryTags['right'], boundaryTags['top'], boundaryTags['bottom'], boundaryTags['front'], boundaryTags['back']]:
        return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        return lambda x,t: fo.hz(x,t)

dirichletConditions = {0:getDBC_hx,
                       1:getDBC_hy,
                       2:getDBC_hz}

fluxBoundaryConditions = {0:'noFlow',
                          1:'noFlow',
                          2:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{},
                                   2:{}}

def stress_u(x,flag):
    if flag not in [boundaryTags['left'],
                    boundaryTags['right']]:
        return 0.0

def stress_v(x,flag):
    if flag not in [boundaryTags['front'],
                    boundaryTags['back']]:
        return 0.0

def stress_w(x,flag):
    if flag not in [boundaryTags['top'],
                    boundaryTags['bottom']]:
        return 0.0

stressFluxBoundaryConditions = {0:stress_u,
                                1:stress_v,
                                2:stress_w}
