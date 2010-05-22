from pyadh import *
from pyadh.default_p import *
from twp_darcy_ls_so_waterTable_2d_p import waterTable,zB,dp,nd

nd = 2
## \todo im_ls_so_darcy_waterTable_2d_p.py doc, push coefficients into TransportCoefficients.py if mature
#

class DarcyNCLevelSetCoefficients(NCLevelSetCoefficients):
    def __init__(self):
        NCLevelSetCoefficients.__init__(self)
    def attachModels(self,modelList):
        self.q_v = modelList[self.flowModelIndex].q[('velocity',0)]
        self.ebq_v = modelList[self.flowModelIndex].ebq[('velocity',0)]
    #for debugging
    def evaluate(self,t,c):
        if c[('dH',0,0)].shape == self.q_v.shape:
            v = self.q_v
        else:
            v = self.ebq_v
        #mwf debug
        #print """Darcy LS v= %s \n""" % v
        self.ncLevelSetCoefficientsEvaluate(v,
                                            c[('u',0)],
                                            c[('grad(u)',0)],
                                            c[('m',0)],
                                            c[('dm',0,0)],
                                            c[('H',0)],
                                            c[('dH',0,0)])
#end DarcyNC

coefficients = DarcyNCLevelSetCoefficients()
coefficients.variableNames=['phi']

class PerturbedWaterTable_phi:
    def __init__(self,waterTable=0.5,zB=0.0,dp=1.0):
        self.zB    =zB
        self.dp    = dp
        self.waterTablePert=lambda x : self.dp*x*(1.0-x) + waterTable + self.zB
    def uOfX(self,X):
        zP = self.waterTablePert(X[0])
        return zP-X[1]
    #end
    def uOfXT(self,X,t):
        return self.uOfX(X)
    #end
#end PerturbedWaterTable

analyticalSolutions = None

def getDBC(x):
    pass

dirichletConditions = {0:getDBC}
initialConditions  = {0:PerturbedWaterTable_phi(waterTable,zB,dp)}

fluxBoundaryConditions = {0:'noFlow'} #for HJ

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

T = 1.0e1
