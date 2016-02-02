from proteus import *
from proteus.default_p import *
from floating_bar import *
from proteus.mprans import Kappa

LevelModelType = Kappa.LevelModel
if useOnlyVF:
    RD_model = None
    LS_model = None
    dissipation_model = 3
    ME_model = 2
else:
    RD_model = 3
    LS_model = 2
    ME_model = 5
    dissipation_model = 6
if movingDomain:
    dissipation_model += 1
    ME_model += 1
#
dissipation_model_flag = 1 
if useRANS == 2:
    dissipation_model_flag=2
elif useRANS == 3:
    dissipation_model_flag=3

coefficients = Kappa.Coefficients(V_model=0,ME_model=ME_model,LS_model=LS_model,RD_model=RD_model,dissipation_model=dissipation_model,
                                  dissipation_model_flag=dissipation_model_flag,#1 -- K-epsilon, 2 -- K-omega 1998, 3 -- K-omega 1988
                                  useMetrics=useMetrics,
                                  rho_0=rho_0,nu_0=nu_0,
                                  rho_1=rho_1,nu_1=nu_1,
                                  g=g,
                                  c_mu=0.09,sigma_k=1.0, 
                                  sc_uref=kappa_sc_uref,
                                  sc_beta=kappa_sc_beta)


def getDBC_k(x,flag):
    if flag == boundaryTags['left']:
        return lambda x,t:kInflow
    if flag == boundaryTags['obstacle']:
        return lambda x,t:0.0
    if flag in [boundaryTags['front'], boundaryTags['back']]:
        if openSides:
            return lambda x,t:kInflow
dirichletConditions = {0:getDBC_k}
#fluxBoundaryConditions = {0:'outFlow'}

def getAFBC_k(x,flag):
    if flag == boundaryTags['left']:
        return None#weak Dirichlet
    elif flag == boundaryTags['obstacle']:
        return None#weak Dirichlet
    elif flag == boundaryTags['right']:
        return None#outflow
    elif flag == boundaryTags['top']:
        if openTop:
            return None#outflow
        else:
            return lambda x,t: 0.0
    elif flag in [boundaryTags['front'], boundaryTags['back']]:
        if openSides:
            return None#outflow
        else:
            return lambda x,t: 0.0
    elif flag == boundaryTags['bottom']:
        return None#outflow
def getDFBC_k(x,flag):
    if flag == boundaryTags['left']:
        return None#weak Dirichlet
    elif flag == boundaryTags['obstacle']:
        return None#weak Dirichlet
    elif flag == boundaryTags['right']:
        return lambda x,t: 0.0 #outflow
    elif flag in [boundaryTags['front'],boundaryTags['back'],boundaryTags['top'],boundaryTags['bottom']]:
        return lambda x,t: 0.0#outflow or no flow


advectiveFluxBoundaryConditions =  {0:getAFBC_k}
diffusiveFluxBoundaryConditions = {0:{0:getDFBC_k}}



class ConstantIC:
    def __init__(self,cval=0.0):
        self.cval=cval
    def uOfXT(self,x,t):
        return self.cval
   
initialConditions  = {0:ConstantIC(cval=kInflow*0.001)}
