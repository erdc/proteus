from pyadh import *
from pyadh.default_p import *
from pyadh.TransportCoefficients import *
from impes_modelParams_forsyth2_2d import *

phase = 'potential' # 'saturation' or 'potential'

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

def getDBC_pressure_HydroRight(x):
    if x[0] == L[0]:
        return lambda x,t: (x[1]-waterTableHeight)*dimensionless_gravity[1]*dimensionless_density

def getDBC_pressure_Box(x):
    pass

def getRecharge_2D_UpperLeft(x):
    if x[1] == L[1]:
        if x[0] < rechargeXboundary:
            #mwf debug
            #print """setting recharge rate at x=%s """ % x
            return lambda x,t: rechargeRate
def getDummyFlux(x):
    pass
    
#getDBC = getDBC_pressure_HydroRight
getDBC = getDBC_pressure_Box

dirichletConditions = {0:getDBC}

class ConstIC_2D_pressure:
    def uOfXT(self,x,t):
        return initialPressure
class LinearIC_2D_pressure:
    def uOfXT(self,x,t):
        return (pondingPressure-0.0)/L[1]*(x[1]-0.0) + 0.0
class HydroIC_2D_pressure:
    def uOfXT(self,x,t):
        return (x[1]-waterTableHeight)*dimensionless_gravity[1]*dimensionless_density

initialConditions  = {0:ConstIC_2D_pressure()}
#initialConditions  = {0:HydroIC_2D_pressure()}

def getAFBC_pres(x):
    if x[0] in [0.0,1.0]:
	return lambda x,t: 0.0
#    return lambda x,t:0.0

fluxBoundaryConditions = {0:'setFlow'}
advectiveFluxBoundaryConditions =  {0:getDummyFlux}
diffusiveFluxBoundaryConditions = {0:{0:getRecharge_2D_UpperLeft}}
