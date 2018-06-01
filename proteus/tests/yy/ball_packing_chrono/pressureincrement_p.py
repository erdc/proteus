from math import *
from proteus import *
from proteus.default_p import *
from rockyRiver import *


#domain = ctx.domain
#nd = ctx.nd
name = "pressureincrement"

from proteus.mprans import PresInc
# coefficients=PresInc.Coefficients(rho_f_min = (1.0-1.0e-8)*rho_1,
#                                  rho_s_min = (1.0-1.0e-8)*rho_s,
#                                  nd = nd,
#                                  modelIndex=PINC_model,
#                                  fluidModelIndex=V_model,
#                                  INTEGRATE_BY_PARTS_DIV_U=False)

LevelModelType = PresInc.LevelModel

class My_coefficients(PresInc.Coefficients):
    def __init__(self):
        PresInc.Coefficients.__init__(self,
                                      rho_f_min = (1.0-1.0e-8)*rho_1,
                                      rho_s_min = (1.0-1.0e-8)*rho_s,
                                      nd = nd,
                                      modelIndex=PINC_model,
                                      fluidModelIndex=V_model,
                                      INTEGRATE_BY_PARTS_DIV_U=False,
                                      )
    def preStep(self, t, firstStep=False):
        """
        Move the current values to values_last to keep cached set of values for bdf1 algorithm
        """
        copyInstructions = {}
        return copyInstructions
    def postStep(self, t, firstStep=False):
        """
        Give the TC object an opportunity to modify itself before the time step.
        """
        #put it after the computation of the pressure 
#         alphaBDF = self.fluidModel.timeIntegration.alpha_bdf
#         self.fluidModel.q[('velocity', 0)][..., 0] -= self.model.q[('grad(u)', 0)][..., 0] / (self.rho_f_min * alphaBDF)
#         self.fluidModel.q[('velocity', 0)][..., 1] -= self.model.q[('grad(u)', 0)][..., 1] / (self.rho_f_min * alphaBDF)
#    

        copyInstructions = {}
        return copyInstructions

coefficients = My_coefficients()



#pressure increment should be zero on any pressure dirichlet boundaries
def getDBC_phi(x,flag):
    if flag in [boundaryTags['top']]: 
        return lambda x,t: 0.0

#the advectiveFlux should be zero on any no-flow  boundaries
def getAdvectiveFlux_qt(x,flag):
    if flag == boundaryTags['top']:
        return None 
    else:
        return lambda x,t: 0.0

def getDiffusiveFlux_phi(x,flag):
    if flag == boundaryTags['top']:
        return None 
    else:
        return lambda x,t: 0.0

class getIBC_phi:
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:getIBC_phi()}

dirichletConditions = {0:getDBC_phi}
advectiveFluxBoundaryConditions = {0:getAdvectiveFlux_qt}
diffusiveFluxBoundaryConditions = {0:{0:getDiffusiveFlux_phi}}
