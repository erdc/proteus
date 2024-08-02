from proteus import *
from proteus.default_n import *
from sw_dam_break_1d_p import *

"""
Straight forward RKDG approximation for the 1D shallow water equations
"""

#pick the polynomial order for the trial (and test spaces)
spaceOrder = 1
assert spaceOrder in [0,1]

#pick the target Courant number based on the polynomial order 
if spaceOrder == 1:
    runCFL=0.1
else:
    runCFL=0.15


#SSP Runge Kutta TimeIntegration with the possible use of limiters (filters)
timeIntegration = SSPRKPIintegration
#time integration order
timeOrder = min(spaceOrder+1,3)
nStagesTime=timeOrder
#time step controller picks the time steps based on the Courant number. This one is needed for
#the Explicit Runge-Kutta discretizations
stepController=Min_dt_RKcontroller
#number of time steps

nDTout = 100

#pick the type of limiter to use
if spaceOrder == 1:
    limiterType = DGlimiterP1Lagrange1d
    #limiterType = DGlimiterP1Lagrange1d_Sw

#Pick the Finite Element Space
if spaceOrder == 1:
    femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis,
                 1:DG_AffineLinearOnSimplexWithNodalBasis}
else:
    femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis,
                 1:DG_AffineP0_OnSimplexWithMonomialBasis}

#Element Quadrature Choice
if spaceOrder == 0:
    elementQuadrature = SimplexLobattoQuadrature(nd,1)
    elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
else:
    elementQuadrature = SimplexGaussQuadrature(nd,2*spaceOrder+1)
    elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)


nn = 31
nLevels = 1

#just used to compute the Eigenvalues for Courant number selection 
subgridError = ShallowWater_CFL(coefficients,nd,g)

numericalFluxType = RusanovNumericalFlux#ShallowWaterHLL_1D#RusanovNumericalFlux#ShallowWater_1D

#
multilevelNonlinearSolver=Newton
#if using explicity SSP Runge Kutta, this 'nonlinear solver skips some of the unneeded computations' 
levelNonlinearSolver=Newton
levelNonlinearSolver =SSPRKNewton
usingSSPRKNewton=True

tolFac = 0.0

nl_atol_res = 1.0e-8

#only save the solution at the 'user time steps' in tnList, which is set using nDTout above
archiveFlag = ArchiveFlags.EVERY_USER_STEP