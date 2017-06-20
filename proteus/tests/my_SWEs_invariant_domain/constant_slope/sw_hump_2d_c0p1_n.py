from proteus import *
from proteus.default_n import *
from sw_hump_2d_p import *

refinement=4
runCFL=0.5
#use_EV_stabilization=True
#use_first_order_flatB_GP_stabilization=True
#use_second_order_flatB_GP_stabilization=True
#use_second_order_NonFlatB_GP_stabilization=True
#use_EV_stabilization = True
#use_second_order_NonFlatB_with_EV_stabilization=True
#timeIntegration_sw2d = "SSP33"
timeIntegration_sw2d = "FE"

#timeIntegration = SSP33
#timeIntegration = BackwardEuler_cfl
#timeIntegration = EdgeBased_ForwardEuler
timeIntegration = SW2DCV.RKEV 
stepController = Min_dt_controller
if timeIntegration_sw2d == "SSP33": #mwf hack
    timeOrder = 3
    nStagesTime = 3
else:
    timeOrder = 1
    nStagesTime = 1

rtol_u[0] = 1.0e-4
rtol_u[1] = 1.0e-4
rtol_u[2] = 1.0e-4
atol_u[0] = 1.0e-4
atol_u[1] = 1.0e-4
atol_u[2] = 1.0e-4
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

multilevelNonlinearSolver  = Newton
levelNonlinearSolver = ExplicitLumpedMassMatrixShallowWaterEquationsSolver
#levelNonlinearSolver = ExplicitConsistentMassMatrixShallowWaterEquationsSolver
#levelNonlinearSolver = Newton

fullNewtonFlag = False #NOTE: False just if the method is explicit
nDTout=100

nnx0=6
nnz=1

nnx = (nnx0-1)*(2**refinement)+1
nny = (nnx-1)/10+1
#nny = 2

he = L[0]/float(nnx-1)
triangleOptions="pAq30Dena%f"  % (0.5*he**2,)

#added flag for using SUPG stabilization based on Berger and Stockstill, 95
try_supg_stabilization = True
subgridError = SW2D.SubgridError(coefficients,nd,lag=True)

massLumping=False

shockCapturing = SW2D.ShockCapturing(coefficients,nd,shockCapturingFactor=0.1,lag=True)

numericalFluxType = SW2DCV.NumericalFlux
#numericalFluxType = DoNothing

tolFac = 0.0

levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest             = 'r-true'

maxLineSearches=0
nl_atol_res = 1.0e-5
nl_rtol_res = 0.0
l_atol_res = 1.0e-7
l_rtol_res = 0.0

matrix = SparseMatrix
multilevelLinearSolver = LU
levelLinearSolver = LU

#conservativeFlux = {0:'pwl'}
