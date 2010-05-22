from pyadh import *
from pyadh.default_p import *
from HypPar_level1_TransportCoefficients import *
from HypPar_SetParams import *
#John Chrispell, Spring 08

phase = "ADR"

analyticalSolution = {0:Ern_and_Proft_Test1()}
coefficients = HypPar_Coefficients(q=q,alpha=alpha,delta=delta,phase=phase,rFuncClass=rFuncClass,setParamsFunc=setParams)

#now define the Dirichlet boundary conditions

def getDBC(x):
    if x[0] == 0.0:
        return lambda x,t: Ern_and_Proft_Test1().uOfXT(x,t)
    if x[0] == 2.0:
        return lambda x,t: Ern_and_Proft_Test1().uOfXT(x,t)
def dummyBC(x):
    pass
dirichletConditions = {0:getDBC} #for Riemann problem {0:dummyBC}


class ParabolicIC:
    def __init__(self,uleft=0.0,uright=0.0):
        self.uLeft=uleft; self.uRight=uright
    def uOfXT(self,x,t):
        if (x[0] >= 0.0 and x[0] <= 1.0):
	     return 16.0*x[0]*x[0]*(1.0-x[0])**2.0 
	else: 
	     return 0.0
class ParabolicIC2:
    def __init__(self,uleft=0.0,uright=0.0):
        self.uLeft=uleft; self.uRight=uright
    def uOfXT(self,x,t):
        if (x[0] >= 0.0 and x[0] <= 1.0):
	     return 16.0*x[0]*x[0]*(1.0-x[0])**2.0 
	else: 
	     return 0.0
	     
class ParabolicIC3:
    def __init__(self,uleft=0.0,uright=0.0):
        self.uLeft=uleft; self.uRight=uright
    def uOfXT(self,x,t):
        if (x[0] >= 0.2 and x[0] <= 1.0):
	     return 1.0
	else: 
	     return 0.0     
	     
	     

#initialConditions  = {0:ShockIC(getDBC)}#{0:RiemIC(uleft=1.0,uright=0.0)}#{0:RiemSlopeIC(uleft=1.0,uright=0.0,du=0.02)}
#initialConditions  = {0:ParabolicIC3(getDBC)}

initialConditions  = {0:Ern_and_Proft_Test1()}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{}}

