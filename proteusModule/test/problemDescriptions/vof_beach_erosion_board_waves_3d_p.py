from pyadh import *
from pyadh.default_p import *
from beach_erosion_board_waves_3d import *
from pyadh import VOF,VolumeAveragedVOF

if useVOF:
    if useSpongeLayer:
        LevelModelType = VOF.OneLevelVOF
    else:
        LevelModelType = VolumeAveragedVOF.OneLevelVolumeAveragedVOF

if useSpongeLayer:
    coefficients = VolumeAveragedVOFCoefficients(LS_model=1,V_model=0,RD_model=3,ME_model=2,epsFact=epsFact_vof,
                                                 setParamsFunc=spongeLayerFunc,checkMass=checkMass)
else:
    coefficients = VOFCoefficients(LS_model=1,V_model=0,RD_model=3,ME_model=2,epsFact=epsFact_vof,checkMass=checkMass)

dirichletConditions = {0:getDBC_vof}

initialConditions  = {0:Flat_H()}

fluxBoundaryConditions = {0:'outFlow'}

advectiveFluxBoundaryConditions =  {0:getAFBC_vof}

diffusiveFluxBoundaryConditions = {0:{}}
