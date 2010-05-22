from pyadh import *
from pyadh.default_p import *
from pyadh.TwophaseDarcyCoefficients import *
from darcy_simp_modelParams import *

phase = 'saturation'

coefficients = TwophaseDarcy_split_saturation(g=g, 
                                              rhon=rhon,
                                              rhow=rhow,
                                              mun    = mun,
                                              muw    = muw,
                                              Ksw=Ksw,
                                              psk_model=model,
                                              omega  = omega,
                                              Sw_max = sw_max,
                                              Sw_min = sw_min,
                                              qScalarConstant = q,
                                              capillaryDiffusionScaling=capillaryDiffusionScaling)

#now define the Dirichlet boundary conditions
riemFudgeFactor = 0.005
uLeft = Se_top*(sw_max-sw_min)+sw_min
uRight= Se_bottom*(sw_max-sw_min)+sw_min
def getDBC_saturation(x,flag):
    if x[0] == 0.0:
        return lambda x,t: uLeft

getDBC = getDBC_saturation

dirichletConditions = {0:getDBC}

class ShockIC:
    def __init__(self,dbc):
        pass
    def uOfXT(self,x,t):
        if x[0] <= riemFudgeFactor:
            return uLeft
	else:
	    return uRight

initialConditions  = {0:ShockIC(getDBC)}
analyticalSolution = {0:AnalyticalSolutions.Buckley_Leverett_RiemannSoln(coefficients,uLeft=uLeft,uRight=uRight,
                                                                         t0=0.0,x0=riemFudgeFactor,T=T,useShallowCopyCoef=False)}
fluxBoundaryConditions = {0:'outFlow'}
advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}


