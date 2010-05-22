from pyadh import *
from pyadh.default_p import *
from pyadh.TransportCoefficients import *
from impes_modelParams_2d import *

phase = 'potential' # 'saturation' or 'potential'

analyticalSolutions = None
coefficients = TwoPhaseFlow(q=q,
                            a=a,
			    Ksw=Ksw,
			    rhon=rhon,
			    rhow=rhow,
			    g=g,
			    mvg_alpha=mvg_alpha,
			    bc_lambda=bc_lambda,
			    bc_pd = bc_pd,
			    mvg_n = mvg_n,
			    mvg_m = mvg_m,
			    omega=omega,
			    mun=mun,
			    muw=muw,
			    sw_max=sw_max,
			    sw_min=sw_min,
			    M=M,
			    R=R,
			    Temp=Temp,
			    p_o=p_o,
			    model=model,
			    phase=phase)

#now define the Dirichlet boundary conditions

def getDBC_potential(x):
    if x[1] == 0.0:
    	return lambda x,t: 0.1
    if x[1] == 1.0:
    	return lambda x,t: 0.0

getDBC = getDBC_potential

dirichletConditions = {0:getDBC}

class DummyIC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0
	
initialConditions  = {0:DummyIC()}

fluxBoundaryConditions = {0:'outFlow'}

def getAFBC_pres(x):
    if x[0] in [0.0,1.0]:
	return lambda x,t: 0.0
#    return lambda x,t: 0.0
#   if x[1] == 0.0:
#	return lambda x,t: -q
    #if x[1] == 1.0:
	#return lambda x,t: -q

advectiveFluxBoundaryConditions =  {0:getAFBC_pres}
diffusiveFluxBoundaryConditions = {0:{}}


