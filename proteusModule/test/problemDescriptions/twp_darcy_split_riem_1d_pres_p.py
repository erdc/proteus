from pyadh import *
from pyadh.default_p import *
from pyadh.TwophaseDarcyCoefficients import *
from twp_darcy_riem_1d import *

phase = 'potential' 

analyticalSolutions = None
coefficients = TwophaseDarcy_split_pressure(g=g, 
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

def getDBC_potential(x,flag):
    if x[0] == 0.0:
    	return lambda x,t: psi_top
    if x[0] == L[0]:
    	return lambda x,t: psi_bottom

getDBC = getDBC_potential

dirichletConditions = {0:getDBC}

class DummyIC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0
	
initialConditions  = {0:DummyIC()}

fluxBoundaryConditions = {0:'outFlow'}

def getAFBC(x):
    pass
#    if x[0] == 0.0:
#	return lambda x,t: -q

advectiveFluxBoundaryConditions =  {0:getAFBC}
diffusiveFluxBoundaryConditions = {0:{}}


