from pyadh import *
from pyadh.default_n import *
from ladr_2d_p import *
timeIntegration = BackwardEuler
femSpaces = {0:DG_AffineP3_OnSimplexWithMonomialBasis}
elementQuadrature = SimplexGaussQuadrature(nd,5)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,5)
nnx=11; nny=11
tnList=[float(i)/10.0 for i in range(11)]
matrix = SparseMatrix
multilevelLinearSolver = LU
timeOrder = 3
nStagesTime = timeOrder
runCFL = 0.125
limiterType = TimeIntegration.DGlimiterPkMonomial2d#

timeIntegration = SSPRKPIintegration
usingSSPRKNewton = True
#end wrapper
#BackwardEuler, ForwardEuler
stepController=Min_dt_RKcontroller
numericalFluxType = Advection_DiagonalUpwind_Diffusion_LDG

archiveFlag = ArchiveFlags.EVERY_USER_STEP
