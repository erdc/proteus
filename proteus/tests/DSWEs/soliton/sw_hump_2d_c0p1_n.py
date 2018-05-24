from proteus import *
from proteus.default_n import *
from sw_hump_2d_p import *

refinement=5
runCFL=0.1
SSPOrder=1

multilevelNonlinearSolver  = Newton
if (LUMPED_MASS_MATRIX==1):
    levelNonlinearSolver = ExplicitLumpedMassMatrixShallowWaterEquationsSolver
else:
    levelNonlinearSolver = ExplicitConsistentMassMatrixShallowWaterEquationsSolver

timeIntegration = DSW2DCV.RKEV
stepController = Min_dt_controller
timeOrder = SSPOrder
nStagesTime = SSPOrder

rtol_u[0] = 1.0e-4
rtol_u[1] = 1.0e-4
rtol_u[2] = 1.0e-4
atol_u[0] = 1.0e-4
atol_u[1] = 1.0e-4
atol_u[2] = 1.0e-4
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis,
             3:C0_AffineLinearOnSimplexWithNodalBasis,
             4:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3) #3
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

fullNewtonFlag = False #NOTE: False just if the method is explicit


nnx0=6
nnz=1

nnx = (nnx0-1)*(2**refinement)+1
nny = (nnx-1)/10+1
#nny = 2

he = L[0]/float(nnx-1)
triangleOptions="pAq30Dena%f"  % (0.5*he**2,)

#added flag for using SUPG stabilization based on Berger and Stockstill, 95
try_supg_stabilization = True
subgridError = DSW2DCV.SubgridError(coefficients,nd,lag=True)

massLumping=False

shockCapturing = DSW2DCV.ShockCapturing(coefficients,nd,shockCapturingFactor=0.1,lag=True)
numericalFluxType = DSW2DCV.NumericalFlux

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
tnList=[0.,1E-6]+[float(n)*T/float(nDTout) for n in range(1,nDTout+1)]
