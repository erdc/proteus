from pyadh import *
from pyadh.default_p import *
from pyadh.TransportCoefficients import *
from impes_modelParams_forsyth2_2d import *

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
			    sw_min = sw_min, 
			    sw_max = sw_max, 
			    model=model,
			    phase=phase,
			    setParamsFunc=setParams)

#now define the Dirichlet boundary conditions

# The boundary conditions here are still not set correctly. 
#  JCC  - 7/3/07


FudgeFactor = 1.0e-1

def getDBC_saturation(x):    
    #constant saturation over slit
    if((x[1] == L[1]) and (x[0] < rechargeXboundary)):
            return lambda x,t: 1.0-FudgeFactor
def getDummyFlux(x):
    pass

def getRecharge_2D_UpperLeft(x):
    if x[1] == L[1]:
        if x[0] < rechargeXboundary:
            #mwf debug
            #print """setting recharge rate at x=%s """ % x
            return lambda x,t: rechargeRate

getDBC = getDBC_saturation

dirichletConditions = {0:getDBC}

class ConstIC_2D_saturation:
    def uOfXT(self,x,t):  
        if ((x[1] == L[1]) and (x[0] < rechargeXboundary)):
            return 1.0 - FudgeFactor
        else:
            return 0.0 + FudgeFactor
    
initialConditions  = {0:ConstIC_2D_saturation()}

#def getAFBC_sat(x):
#    if x[0] in [0.0,1.0]:
#	return lambda x,t: 0.0
#    return lambda x,t:0.0

fluxBoundaryConditions = {0:'noFlow'}
advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}



