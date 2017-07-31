from proteus.default_n import *
from proteus import (Context,)
from proteus.mprans import SW2D
import sw_hump_2d_p
Context.setFromModule(sw_hump_2d_p)
ct = Context.get()

reflecting_BCs=ct.opts.reflecting_BCs
refinement=ct.opts.refinement
runCFL=0.5
timeIntegration_sw2d = "SSP33"
#timeIntegration_sw2d = "FE"

multilevelNonlinearSolver  = Newton
if (ct.LUMPED_MASS_MATRIX==1):
    levelNonlinearSolver = ExplicitLumpedMassMatrixShallowWaterEquationsSolver
else:
    levelNonlinearSolver = ExplicitConsistentMassMatrixShallowWaterEquationsSolver


timeIntegration = ct.SW2DCV.RKEV 
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

elementQuadrature = SimplexGaussQuadrature(ct.nd,3)
elementBoundaryQuadrature = SimplexGaussQuadrature(ct.nd-1,3)

fullNewtonFlag = False #NOTE: False just if the method is explicit


nnx0=6
nnz=1

nnx = (nnx0-1)*(2**refinement)+1
nny = (nnx-1)/2+1
#nny = 2

he = ct.L[0]/float(nnx-1)
triangleOptions="pAq30Dena%f"  % (0.5*he**2,)

#added flag for using SUPG stabilization based on Berger and Stockstill, 95
try_supg_stabilization = True
subgridError = SW2D.SubgridError(ct.coefficients,ct.nd,lag=True)

massLumping=False

shockCapturing = SW2D.ShockCapturing(ct.coefficients,ct.nd,shockCapturingFactor=0.1,lag=True)

numericalFluxType = SW2D.NumericalFlux

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
tnList=[0.,1E-6]+[float(n)*ct.T/float(ct.nDTout) for n in range(1,ct.nDTout+1)]

