from pyadh import *
from pyadh.default_p import *
from threep_cylinder_md_2d import *

initialConditions = None

analyticalSolution = {}

coefficients = MovingMesh(E=10.0,nu=0.3,g=[0.0,0.0],nd=nd,moveMesh=movingDomain)

fixedBoundary =  [boundaryTags['bottom'],boundaryTags['top'],boundaryTags['upstream'],boundaryTags['downstream']]

class TranslatingObstacle(AuxiliaryVariables.AV_base):
    def __init__(self):
        self.current_center=cylinder_center
        self.r = 0.5*cylinder_radius
        self.omega=1.0/1.0
        self.cx = cylinder_center[0]+self.r
        self.cy = cylinder_center[1]
    def attachModel(self,model,ar):
        self.model=model
        return self
    def center(self,t):
        return (self.cx+self.r*cos(self.omega*2.0*math.pi*t+math.pi),
                self.cy+self.r*sin(self.omega*2.0*math.pi*t +math.pi))
    def hx(self,t):
        h = tro.center(t)[0]-tro.current_center[0]
        return h
    def hy(self,t):
        h = tro.center(t)[1]-tro.current_center[1]
        return h
    def calculate(self):
        self.current_center=self.center(self.model.stepController.t_model_last)

class TranslatingFreeObstacle(AuxiliaryVariables.AV_base):
    def __init__(self):
        self.current_center=navier_stokes_cylinder_2d_p.cylinder_center
        self.r = 0.5*navier_stokes_cylinder_2d_p.cylinder_radius
        self.omega=1.0/1.0
        self.cx = navier_stokes_cylinder_2d_p.cylinder_center[0]+self.r
        self.cy = navier_stokes_cylinder_2d_p.cylinder_center[1]
        self.body=None
        self.h=(0.0,0.0)
    def attachModel(self,model,ar):
        self.model=model
        return self
    def attachAuxiliaryVariables(self,avDict):
        self.object = avDict.values()[0][0]
    def hx(self,t):
        return self.object.h[0]
        #return self.h[0]
    def hy(self,t):
        return self.object.h[1]
        #return self.h[1]
    def calculate(self):
        self.h=(self.object.position[0]-self.object.last_position[0],
                self.object.position[1]-self.object.last_position[1])
        self.h=self.object.h#(self.object.position[0]-self.object.last_position[0],
                #self.object.position[1]-self.object.last_position[1])
        print "displacement------------------------------------",self.h[0],self.h[1]
tro = TranslatingFreeObstacle()

def getDBC_u(x,flag):
    if flag in fixedBoundary:
        return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        return lambda x,t: tro.hx(t)

def getDBC_v(x,flag):
    if flag in fixedBoundary:
        return lambda x,t: 0.0
    if flag == boundaryTags['obstacle']:
        return lambda x,t: tro.hy(t)

dirichletConditions = {0:getDBC_u,
                       1:getDBC_v}

fluxBoundaryConditions = {0:'noFlow',
                          1:'noFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{},1:{}}

