from pyadh import *
from pyadh.default_n import *
from re_minesCell1_2d_p import *

#timeIntegrator  = ForwardIntegrator
#timeIntegration = BackwardEuler
timeIntegration = FLCBDF
stepController = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
rtol_u[0] = 1.0e-2
atol_u[0] = 1.0e-2
runCFL=None
DT = None
nDTout = 1000

triangleOptions= "q30DenA" #default "q30Den"
#triangleOptions+="A"# for element type data
nn=3
nLevels = 4

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}

#appears to work
elementQuadrature = {'default':SimplexLobattoQuadrature(nd,3),
                     ('u',0):SimplexGaussQuadrature(nd,3)}
elementBoundaryQuadrature ={'default':SimplexGaussQuadrature(nd-1,2),
                            ('u',0):SimplexGaussQuadrature(nd-1,2)}

#mwf looks like works too
#elementQuadrature = {'default':SimplexGaussQuadrature(nd,3),
#                     ('u',0):SimplexGaussQuadrature(nd,3)}
#elementBoundaryQuadrature = {'default':SimplexGaussQuadrature(nd-1,3),
#                             ('u',0):SimplexGaussQuadrature(nd-1,3)}



subgridError = None

shockCapturing = None

massLumping = False

numericalFluxType = None

#multilevelNonlinearSolver  = NLNI
multilevelNonlinearSolver = Newton

levelNonlinearSolver = Newton

nonlinearSmoother = NLJacobi

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-6

maxNonlinearIts = 10#0#1

maxLineSearches = 10#0
matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = StarILU

linTolFac = 0.001

conservativeFlux = {0:'pwl'}#{0:'pwc'}#{0:'sun-gs-rt0'}#
#if using sun-gs-rt0 need weak dirichlet conditions
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior

