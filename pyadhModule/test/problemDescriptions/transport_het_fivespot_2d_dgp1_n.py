from pyadh import *
from pyadh.default_n import *
from transport_het_fivespot_2d_p import *

timeOrder = 2
nStagesTime = timeOrder
runCFL = 0.15
timeIntegration = SSPRKPIintegration#LinearSSPRKintegration
limiterType = DGlimiterP1Lagrange2d
stepController = Min_dt_RKcontroller
#timeOrder = 1
#timeIntegration = BackwardEuler_cfl
#stepController = Min_dt_controller

DT = None
femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,gw_quad_order)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,gw_quad_order)

nLevels = 1
massLumping=False
#subgridError = None
subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton
fullNewtonFlag = True

tolFac = 1.0e-8

nl_atol_res = 1.0e-8

matrix = SparseMatrix
#matrix = Numeric.array

multilevelLinearSolver = LU
#multilevelLinearSolver = NI
numericalFluxType = Advection_DiagonalUpwind
#Advection_DiagonalUpwind_Diffusion_IIPG

levelLinearSolver = LU
#levelLinearSolver = MGM
#levelLinearSolver = StarILU
#levelLinearSolver = GaussSeidel
#levelLinearSolver = Jacobi

smoother = StarILU
smoother = GaussSeidel
smoother = Jacobi

linTolFac = 1.0e-10

conservativeFlux = None

archiveFlag = ArchiveFlags.EVERY_USER_STEP
