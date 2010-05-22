from pyadh import *
from pyadh.default_p import *
from pyadh.TransportCoefficients import *
from impes_modelParams_finger_2d import *

phase = 'saturation'

analyticalSolutions = None
coefficients = TwoPhaseFlowHet(q=q,
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
			    sw_min=sw_min,
			    sw_max=sw_max,
			    model=model,
			    phase=phase,
			    setParamsFunc=setParams)
			    
    

#now define the Dirichlet boundary conditions

def getDBC_saturation(x):
    if x[1] == 0.0:
       return lambda x,t: 1.0 - FudgeFactor
    if x[1] == 1.0:
        return lambda x,t: FudgeFactor
getDBC = getDBC_saturation

dirichletConditions = {0:getDBC}

class ShockIC:
    def __init__(self,dbc):
        pass
    def uOfXT(self,x,t):
        Amp = 0.2
	x_c = 0.5
	dx  = 0.2
	shockLine = 0.4
	if x[1] <= shockLine:
	    return 1.0 - FudgeFactor
        elif ((x[1] <= shockLine + Amp*cos( math.pi*(x[0]-x_c)/dx) ) and ( abs(x[0] - x_c) < dx/2 )):
	    return 1.0 - FudgeFactor
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


