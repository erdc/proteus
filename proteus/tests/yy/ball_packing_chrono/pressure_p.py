from math import *
from proteus import *
from proteus.default_p import *
from rockyRiver import *
from proteus.mprans import Pres

name = "pressure"
LevelModelType = Pres.LevelModel

# coefficients=Pres.Coefficients(modelIndex=PRESSURE_model,
#                                fluidModelIndex=V_model,
#                                pressureIncrementModelIndex=PINC_model,
#                                useRotationalForm=False)
class My_coefficients(Pres.Coefficients):
    def __init__(self):
        Pres.Coefficients.__init__(self,
                                  modelIndex=PRESSURE_model,
                                  fluidModelIndex=V_model,
                                  pressureIncrementModelIndex=PINC_model,
                                  useRotationalForm=True)
    
    def attachModels(self, modelList):
        Pres.Coefficients.attachModels(self, modelList)
        
        self.q_old         = np.copy(self.model.q[('u', 0)])
        self.ebqe_old      = np.copy(self.model.ebqe[('u', 0)])
        self.q_grad_old    = np.copy(self.model.q[('grad(u)', 0)])
        self.ebqe_grad_old = np.copy(self.model.ebqe[('grad(u)', 0)])
        
    def postStep(self, t, firstStep=False):
        
        if firstStep==True:##since this is the first step and prestep do nothing
            self.q_old[:] = self.model.q[('u', 0)]
            self.ebqe_old[:] = self.model.ebqe[('u', 0)]
            self.q_grad_old[:] = self.model.q[('grad(u)', 0)]
            self.ebqe_grad_old[:] = self.model.ebqe[('grad(u)', 0)]
            
        
        self.model.q_p_sharp[:] = self.model.q[('u', 0)]
        self.model.ebqe_p_sharp[:] = self.model.ebqe[('u', 0)]
        self.model.q_grad_p_sharp[:] = self.model.q[('grad(u)', 0)]
        self.model.ebqe_grad_p_sharp[:] = self.model.ebqe[('grad(u)', 0)]

        # prediction pressure: 
        # this is not right even for fixed time step. Use the above low order method on p_sharp
        # check this implementation with convergence test.
#         self.model.q_p_sharp[:] = 2*self.model.q[('u', 0)]-self.q_old
#         self.model.ebqe_p_sharp[:] = 2*self.model.ebqe[('u', 0)]-self.ebqe_old
#         self.model.q_grad_p_sharp[:] = 2*self.model.q[('grad(u)', 0)]-self.q_grad_old
#         self.model.ebqe_grad_p_sharp[:] = 2*self.model.ebqe[('grad(u)', 0)]-self.ebqe_grad_old

        self.q_old[:] = self.model.q[('u', 0)]
        self.ebqe_old[:] = self.model.ebqe[('u', 0)]
        self.q_grad_old[:] = self.model.q[('grad(u)', 0)]
        self.ebqe_grad_old[:] = self.model.ebqe[('grad(u)', 0)]

        ##update velocity since for rotation form the tilde(u) should be used
        alphaBDF = self.fluidModel.timeIntegration.alpha_bdf
        self.fluidModel.q[('velocity', 0)][..., 0] -= self.pressureIncrementModel.q[('grad(u)', 0)][..., 0] / (self.pressureIncrementModel.coefficients.rho_f_min * alphaBDF)
        self.fluidModel.q[('velocity', 0)][..., 1] -= self.pressureIncrementModel.q[('grad(u)', 0)][..., 1] / (self.pressureIncrementModel.coefficients.rho_f_min * alphaBDF)
        
        copyInstructions = {}
        return copyInstructions

coefficients = My_coefficients()
#===============================================================================
# BC
#===============================================================================
def getDBC_p(x,flag):
    if flag in [boundaryTags['top']]:
        return lambda x,t: 0.0

def getFlux(x,flag):
    if not(flag == boundaryTags['top']):
        return lambda x,t: 0.0

class getIBC_p:
    def uOfXT(self,x,t):
        return -(L[1] - x[1])*rho_1*g[1]

initialConditions = {0:getIBC_p()}

dirichletConditions = {0:getDBC_p } # pressure bc are explicitly set
advectiveFluxBoundaryConditions = {0:getFlux}
