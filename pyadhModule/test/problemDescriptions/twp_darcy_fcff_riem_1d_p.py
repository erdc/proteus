from pyadh import *
from pyadh.default_p import *
from pyadh.TwophaseDarcyCoefficients import *
from twp_darcy_riem_1d import *

coefficients = TwophaseDarcy_fc_ff(g=g, 
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

def getDBC_potential(x,flag):
    if x[0] == 0.0:
    	return lambda x,t: psi_top
    if x[0] == L[0]:
    	return lambda x,t: psi_bottom


dirichletConditions = {0:getDBC_saturation,1:getDBC_potential}

class ShockIC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        if x[0] <= 0.5:
            return Se_top*(sw_max-sw_min)+sw_min
	else:
	    return Se_bottom*(sw_max-sw_min)+sw_min

class ZeroIC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions  = {0:ShockIC(),1:ZeroIC()}

fluxBoundaryConditions = {0:'outFlow',1:'outFlow'}
advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}


