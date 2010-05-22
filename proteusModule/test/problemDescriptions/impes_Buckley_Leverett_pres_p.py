from pyadh import *
from pyadh.default_p import *
from pyadh.TwophaseDarcyCoefficients import *
from darcy_simp_modelParams import *

phase = 'potential' 

coefficients = TwophaseDarcy_split_pressure(g=g, 
                                            rhon=rhon,
                                            rhow=rhow,
                                            mun    = mun,
                                            muw    = muw,
                                            Ksw=Ksw,
                                            psk_model=model,
                                            omega  = omega,
                                            Sw_max = sw_max,
                                            Sw_min = sw_min,
                                            capillaryDiffusionScaling=capillaryDiffusionScaling)


#now define the Dirichlet boundary conditions

def getDBC_potential(x,flag):
    #if x[0] == 0.0:
    #	return lambda x,t: 0.1
    if x[0] == L[0]:
    	return lambda x,t: 0.0

getDBC = getDBC_potential

dirichletConditions = {0:getDBC}

class DummyIC:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0
	
initialConditions  = {0:DummyIC()}

fluxBoundaryConditions = {0:'setFlow'}

def getAFBC(x,tag):
    #pass
    if x[0] == 0.0:
	return lambda x,t: -q

advectiveFluxBoundaryConditions =  {0:getAFBC}
diffusiveFluxBoundaryConditions = {0:{}}


