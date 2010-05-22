from pyadh import *
from pyadh.default_p import *
from pyadh.TransportCoefficients import *
from impes_modelParams_2d import *

phase = 'saturation'

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

def getDBC_saturation(x):
    if x[1] == 0.0:
       return lambda x,t: 1.0 - FudgeFactor
    if x[1] == 1.0:
        return lambda x,t: 0.0 + FudgeFactor
getDBC = getDBC_saturation

dirichletConditions = {0:getDBC}

class ShockIC:
    def __init__(self,dbc):
        pass
    def uOfXT(self,x,t):
        if x[1] <= 0.5: #+self.a*sin(2.0*math.pi*x[0]*self.k):
            return 1.0-FudgeFactor
	else:
	    return FudgeFactor

initialConditions  = {0:ShockIC(getDBC)}

def getAFBC_sat(x):
    if x[0] in [0.0,1.0]:
	return lambda x,t: 0.0
#    return lambda x,t:0.0

fluxBoundaryConditions = {0:'outFlow'}
advectiveFluxBoundaryConditions =  {0:getAFBC_sat}
diffusiveFluxBoundaryConditions = {0:{}}


