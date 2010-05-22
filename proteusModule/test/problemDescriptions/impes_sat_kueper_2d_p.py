from pyadh import *
from pyadh.default_p import *
from pyadh.TransportCoefficients import *
from impes_kueper_setParams import *

phase = 'saturation'

analyticalSolutions = None

if useHet:
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
else:
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
				sw_max = sw_max,
				sw_min = sw_min, 
                                model=model,
                                phase=phase)

#now define the Dirichlet boundary conditions

def getDBC_saturation(x):
    #constant saturation over slit
    if x[1] == top:
        if (slit_is_top or
            (x[0] >= right/3.0 and
             x[0] <= 2.0*right/3.0)):
            return lambda x,t: 1.0-FudgeFactor
    if open_bottom:
        if x[1] == 0.0:
            return lambda x,t: 1.0-FudgeFactor
    else:#open sides
        if x[0] in [0.0,right]:
            return lambda x,t: FudgeFactor

getDBC = getDBC_saturation

dirichletConditions = {0:getDBC}

class ShockIC:
    def uOfXT(self,x,t): 
        if (x[1] == top and
            (slit_is_top or
            (x[0] >= right/3.0 and
             x[0] <= 2.0*right/3.0))):
            return 1.0-FudgeFactor
        else:
            return FudgeFactor

initialConditions  = {0:ShockIC()}

fluxBoundaryConditions = {0:'noFlow'}

# def getAFBC(x):
#     #no flow outside of slit and on bottom
#     if x[1] == top:
#         if (x[0] < right/3.0 or
#             x[0] > 2.0*right/3.0):
#             return lambda x,t: 0.0
#     if x[1] == 0.0:
#         return lambda x,t: 0.0
# advectiveFluxBoundaryConditions =  {0:getAFBC}
advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}


