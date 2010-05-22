from pyadh import *
from pyadh.default_p import *
from Lin_Liu_waves import *

if useSpongeLayer == True:
    coefficients = VolumeAveragedVOFCoefficients(LS_model=1,V_model=0,RD_model=3,ME_model=2,
                                                 setParamsFunc=spongeLayerFunc)
else:
    coefficients = VOFCoefficients(LS_model=1,V_model=0,RD_model=3,ME_model=2)
#coefficients = VOFCoefficientsWithWaveMaker(LS_model=1,V_model=0,RD_model=3,ME_model=2,
#                                            waveFlag=waveFlag,
#                                            waveHeight=waveHeight,
#                                            waveCelerity=waveCelerity,
#                                            waveFrequency=waveFrequency,
#                                            waveNumber=waveNumber,
#                                            waterDepth=waterLevelBase,
#                                            Omega_s=Omega_s,
#                                            epsFact_source=epsFact_source)

analyticalSolutions = None

dirichletConditions = {0:getDBC_vof}

initialConditions  = {0:Flat_H()}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_vof}

diffusiveFluxBoundaryConditions = {0:{}}
