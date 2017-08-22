from math import *
from proteus import *
from proteus.default_p import *
from dambreak import *
from proteus.mprans import PresInit

#domain = ctx.domain
#nd = ctx.nd
name = "pressureInitial"

coefficients=PresInit.Coefficients(nd=nd,
                                   modelIndex=PINIT_model,
                                   fluidModelIndex=V_model,
                                   pressureModelIndex=PRESSURE_model)

#pressure increment should be zero on any pressure dirichlet boundaries
def getDBC_pInit(x,flag):
    if flag == boundaryTags['top']:
        return lambda x,t: 0.0

#the advectiveFlux should be zero on any no-flow  boundaries
def getAdvectiveFlux_pInit(x,flag):
    if flag != boundaryTags['top']:
        return lambda x,t: 0.0

def getDiffusiveFlux_pInit(x,flag):
    if flag != boundaryTags['top']:
        return lambda x,t: 0.0

class getIBC_pInit:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:getIBC_pInit()}

dirichletConditions = {0:getDBC_pInit }
advectiveFluxBoundaryConditions = {0:getAdvectiveFlux_pInit}
diffusiveFluxBoundaryConditions = {0:{0:getDiffusiveFlux_pInit}}
