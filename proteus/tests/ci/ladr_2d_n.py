from __future__ import absolute_import
from __future__ import division
from builtins import range
from past.utils import old_div
from proteus import *
from proteus.default_n import *
try:
    from .ladr_2d_p import *
except:
    from ladr_2d_p import *

timeIntegration = BackwardEuler_cfl
stepController = Min_dt_cfl_controller
runCFL=1.0
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
elementQuadrature = SimplexGaussQuadrature(nd,3)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,lag=False)
shockCapturing = ResGradQuad_SC(coefficients,nd,
                               shockCapturingFactor=0.99,
                               lag=True)
numericalFluxType = Advection_DiagonalUpwind_Diffusion_SIPG_exterior
nnx=41; nny=41
tnList=[old_div(float(i),40.0) for i in range(11)]
matrix = SparseMatrix
multilevelLinearSolver = LU
linearSmoother = None
l_atol_res = 1.0e-8
parallelPartitioningType = MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 1
conservativeFlux =  None
cfluxtag = None
