from pyadh import *
from pyadh.default_p import *
from pyadh.TransportCoefficients import *
from impes_modelParams_2d import *

analyticalSolutions = None
coefficients = TwophaseDarcyFC(Ksw=Ksw,
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
                               model=model)

#now define the Dirichlet boundary conditions

def getDBC_sw(x):
    if x[1] == 0.0:
        return lambda x,t: FudgeFactor
    if x[1] == L[1]:
        return lambda x,t: 1.0 - FudgeFactor

def getDBC_psiw(x):
    if x[1] == 0.0:
        return lambda x,t: 0.0
    if x[1] == L[1]:
    	return lambda x,t: 0.1

dirichletConditions = {0:getDBC_sw,1:getDBC_psiw}

class sw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if x[1] <= 0.1:
            return FudgeFactor
 	else:
 	    return 1.0-FudgeFactor
#        return (1.0-FudgeFactor)*(L[0]-x[0])+FudgeFactor*x[0]

class psiw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0
	
initialConditions  = {0:sw_IC(),1:psiw_IC()}

fluxBoundaryConditions = {0:'outFlow',1:'outFlow'}


def get_w_AFBC(x):
    if x[0] in [0.0,L[0]]:
        return lambda x,t: 0.0
#    if x[0] == 0.0:
#	return lambda x,t: -q
#    if x[0] == L[0]:
#	return lambda x,t: 0.0

def get_n_AFBC(x):
    if x[0] in [0.0,L[0]]:
        return lambda x,t: 0.0
#    if x[0] == 0.0:
#	return lambda x,t: 0.0
#    if x[0] == L[0]:
#        return lambda x,t: 0.0

advectiveFluxBoundaryConditions =  {0:get_w_AFBC,
                                    1:get_n_AFBC}
diffusiveFluxBoundaryConditions = {0:{},1:{}}


