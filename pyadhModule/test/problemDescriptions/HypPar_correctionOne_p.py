from pyadh import *
from pyadh.default_p import *
from HypPar_TransportCoefficients import *
from HypPar_SetParams import *
#John Chrispell, Spring 08

phase = "correctionOne"

coefficients = HypPar_Coefficients(q=q,alpha=alpha,delta=delta,phase=phase,rFuncClass=rFuncClass,setParamsFunc=setParams)

#now define the Dirichlet boundary conditions

def getDBC(x):
    if x[0] == 0.0:
        return lambda x,t: Ern_and_Proft_Test1().uOfXT(x,t) 
    if x[0] == 4.0:
        return lambda x,t: Ern_and_Proft_Test1().uOfXT(x,t) 
def dummyBC(x):
    pass
dirichletConditions = {0:getDBC}#for Riemann problem {0:dummyBC}
#dirichletConditions = {0:dummyBC}


class ParabolicIC2:
    def __init__(self,uleft=0.0,uright=0.0):
        self.uLeft=uleft; self.uRight=uright
    def uOfXT(self,x,t):
        if (x[0] >= 0.0 and x[0] <= 1.0):
	     return 0.0 
	else: 
	     return 0.0

#initialConditions  = {0:ShockIC(getDBC)}#{0:RiemIC(uleft=1.0,uright=0.0)}#{0:RiemSlopeIC(uleft=1.0,uright=0.0,du=0.02)}
initialConditions  = {0:ParabolicIC2(getDBC)}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

