from proteus import *
from proteus.default_p import *

"""
flow equation for simple miscible displacement example
"""
from miscible_displacement import *

coefficients = MiscibleDisplacementCoefficients_Flow(permeabilities,sources,
                                                     nd=nd,
                                                     viscosity_a=mu_a,
                                                     viscosity_b=mu_b,
                                                     materialValuesLocallyConstant=True)


initialConditions = initialConditions_flow

dirichletConditions = {0:pressure_bc}
fluxBoundaryConditions = {0:'setFlow'}

advectiveFluxBoundaryConditions =  {0:noflux}

diffusiveFluxBoundaryConditions = {0:{0:noflux}}
