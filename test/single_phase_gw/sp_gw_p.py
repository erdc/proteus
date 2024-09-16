import proteus
from proteus import *
from proteus.default_p import *

reload(proteus.default_p)

"""
flow equation for transient, single phase flow example
"""
from .single_phase_gw import *

coefficients = STC.SinglePhaseDarcyCoefficients(conductivities,sources,
                                                S_ss,
                                                nc=1,
                                                nd=nd,
                                                materialValuesLocallyConstant=True,
                                                timeVaryingCoefficients = True)


initialConditions = initialConditions_flow

dirichletConditions = {0:head_bc}
fluxBoundaryConditions = {0:'setFlow'}

advectiveFluxBoundaryConditions =  {0:noflux}

diffusiveFluxBoundaryConditions = {0:{0:noflux}}
