from pyadh import *
from pyadh.default_p import *
from pyadh.TransportCoefficients import *
from impes_kueper_setParams import *

analyticalSolutions = None

if useHet:
    coefficients = TwophaseFFDarcyFCHet(Ksw=Ksw,
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
                                        model=model,
                                        setParamsFunc=setParams)
else:
    coefficients = TwophaseFFDarcyFC(Ksw=Ksw,
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

def getDBC_psiw(x):
    #constant head over slit
    if x[1] == top:
        if (slit_is_top or
            (x[0] >=right/3.0 and
            x[0] <= 2.0*right/3.0)):
            return lambda x,t: 0.1
    if open_bottom:
        if x[1] == 0.0:
            return lambda x,t: 0.0
    else:#open sides
        if x[0] in [0.0,right]:
            return lambda x,t: -x[1]#0.0

dirichletConditions = {0:getDBC_sw,1:getDBC_psiw}

class sw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t): 
        if (x[1] == top and
            (slit_is_top or
            (x[0] >= right/3.0 and
             x[0] <= 2.0*right/3.0))):
            return 1.0-FudgeFactor
        else:
            return FudgeFactor

class psiw_IC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if (x[1] == top and
            (slit_is_top or
            (x[0] >=right/3.0 and
             x[0] <= 2.0*right/3.0))):
            return 0.1
        else:
            return 0.0
	
initialConditions  = {0:sw_IC(),1:psiw_IC()}

fluxBoundaryConditions = {0:'noFlow',1:'noFlow'}


# def get_w_AFBC(x):
#     #no flow outside  of slit and on bottom
#     if x[1] == top:
#         if (not slit_is_top or
#             (x[0] < right/3.0 or
#             x[0] > 2.0*right/3.0)):
#             return lambda x,t: 0.0
#     if x[1] == 0.0:
#         return lambda x,t: 0.0

# def get_n_AFBC(x):
#     #no flow outside  of slit and on bottom
#     if x[1] == top:
#         if (not slit_is_top or
#             (x[0] < right/3.0 or
#             x[0] > 2.0*right/3.0)):
#             return lambda x,t: 0.0
#     if x[1] == 0.0:
#         return lambda x,t: 0.0

# advectiveFluxBoundaryConditions =  {0:get_w_AFBC,
#                                     1:get_n_AFBC}
advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{},1:{}}


