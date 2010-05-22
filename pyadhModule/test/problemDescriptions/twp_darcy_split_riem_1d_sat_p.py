from pyadh import *
from pyadh.default_p import *
from pyadh.TwophaseDarcyCoefficients import *
from twp_darcy_riem_1d import *

phase = 'saturation'

analyticalSolutions = None
coefficients = TwophaseDarcy_split_saturation(g=g, 
                                            rhon=rhon,
                                            rhow=rhow,
                                            mun    = mun,
                                            muw    = muw,
                                            Ksw=Ksw,
                                            psk_model=model,
                                            vg_alpha = mvg_alpha,
                                            vg_m  = mvg_m,
                                            bc_pd  = bc_pd, 
                                            bc_lambda = bc_lambda,
                                            omega  = omega,
                                            Sw_max = sw_max,
                                            Sw_min = sw_min)

#now define the Dirichlet boundary conditions

def getDBC_saturation(x,flag):
    if x[0] == 0.0:
        return lambda x,t: Se_top*(sw_max-sw_min)+sw_min
    if x[0] == L[0]:
        return lambda x,t: Se_bottom*(sw_max-sw_min)+sw_min

getDBC = getDBC_saturation

dirichletConditions = {0:getDBC}

class ShockIC:
    def __init__(self,dbc):
        pass
    def uOfXT(self,x,t):
        if x[0] <= 0.5:
            return Se_top*(sw_max-sw_min)+sw_min
	else:
	    return Se_bottom*(sw_max-sw_min)+sw_min

initialConditions  = {0:ShockIC(getDBC)}

fluxBoundaryConditions = {0:'outFlow'}
advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}


