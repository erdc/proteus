from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from wigley import *
from proteus.mprans import VOF

LevelModelType = VOF.LevelModel
if useOnlyVF:
    RD_model = None
    LS_model = None
else:
    RD_model = 3
    LS_model = 2

coefficients = VOF.Coefficients(LS_model=LS_model,V_model=0,RD_model=RD_model,ME_model=1,
                                checkMass=False,useMetrics=useMetrics,
                                epsFact=epsFact_vof,sc_uref=vof_sc_uref,sc_beta=vof_sc_beta,movingDomain=movingDomain)

def getDBC_vof(x,flag):
    if openTop and flag == boundaryTags['top']:
        return outflowVF
    elif openSides and (flag == boundaryTags['front'] or flag == boundaryTags['back']):
        return outflowVF
    elif flag == boundaryTags['right']:
        if openEnd:
            return outflowVF
        else:
            return None
    elif openSides and (flag == boundaryTags['front'] or flag == boundaryTags['back']):
        return outflowVF
    elif flag == boundaryTags['left']:
        return waveVF

dirichletConditions = {0:getDBC_vof}

def getAFBC_vof(x,flag):
    if flag == boundaryTags['top']:
    	if openTop:
            return None
    	else:
            return lambda x,t: 0.0
    elif (flag == boundaryTags['front'] or flag == boundaryTags['back']):
        if openSides:
            return None
        else:
            return lambda x,t: 0.0
    elif flag == boundaryTags['right']:
        if openEnd:
            return None
        else:
            return lambda x,t: 0.0
    elif flag == boundaryTags['left']:
        return None
    elif openSides and (flag == boundaryTags['front'] or flag == boundaryTags['back']):
        return None
    else:
    	return lambda x,t: 0.0

advectiveFluxBoundaryConditions = {0:getAFBC_vof}
diffusiveFluxBoundaryConditions = {0:{}}

class VF_IC:
    def uOfXT(self,x,t):
        return waveVF_init(x,t)

initialConditions  = {0:VF_IC()}
