from pyadh import *
from pyadh.default_n import *
from ladr_2d_p import *
timeIntegration = FLCBDF
stepController = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
rtol_u[0] = 1.0e-2
atol_u[0] = 1.0e-2

femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=False,stabFlag='1')
shockCapturing=ResGradQuad_SC(coefficients,nd,shockCapturingFactor=0.35,lag=True)
elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)
nnx=11; nny=11
tnList=[float(i)/10.0 for i in range(11)]
matrix = SparseMatrix
multilevelLinearSolver = LU
