from pyadh import *
from pyadh.default_n import *
from poisson_1c_het_2d_p import *

timeIntegration = NoIntegration
nDTout = 1


femSpaces = {0:DG_Constants}

elementQuadrature = {'default':SimplexGaussQuadrature(nd,3)}
elementBoundaryQuadrature = {'default':SimplexGaussQuadrature(nd-1,3)}

nn = 5
nLevels = 3

subgridError = None
#subgridError = AdvectionDiffusionReaction_ASGS(coefficients,nd,stabFlag='2',lag=True)

shockCapturing = None
#shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5e2,lag=True)
#shockCapturing = ResGradJuanes_SC(coefficients,nd,shockCapturingFactor=0.5,lag=False)

numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8

matrix = SparseMatrix
#matrix = Numeric.array

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001

conservativeFlux = {0:'point-eval'}#{0:'sun-gs-rt0'}#{0:'sun-rt0'}#{0:'pwc'}#{0:'pwl'} #{0:'pwl-bdm'} 

archiveFlag = ArchiveFlags.EVERY_USER_STEP
