from pyadh import *
from pyadh.default_n import *
from impes_Buckley_Leverett_sat_p import *

#general type of integration (Forward or to SteadyState)
timeIntegrator  = ForwardIntegrator

limiterType =TimeIntegration.DGlimiterP1Lagrange1d
timeOrder = 2
nStagesTime= timeOrder
timeIntegration = SSPRKPIintegration
usingSSPRKNewton = True
stepController=Min_dt_RKcontroller
runCFL = 0.3




DT=None
nDTout = 1 #int(T/DT)
print "nDTout",nDTout

femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}

#femSpaces = {0:DG_Constants}

elementQuadrature = SimplexGaussQuadrature(nd,3)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn=26
nLevels=1
subgridError = None
massLumping=False
shockCapturing = None

#numericalFluxType = Advection_DiagonalUpwind_Diffusion_LDG
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG
#numericalFluxType = RusanovNumericalFlux_Diagonal_Diffusion_IIPG
numericalFluxType = RusanovNumericalFlux_Diagonal

multilevelNonlinearSolver  = Newton

#levelNonlinearSolver = NLStarILU
levelNonlinearSolver = SSPRKNewton #Newton
maxNonlinearIts = 10

fullNewtonFlag = True

tolFac = 0.0#0.001

nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.0001

conservativeFlux = None

archiveFlag = ArchiveFlags.EVERY_USER_STEP
