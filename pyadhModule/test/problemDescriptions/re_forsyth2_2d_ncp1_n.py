from pyadh import *
from pyadh.default_n import *
from re_forsyth2_2d_p import *

nn=3
#nLevels = 3
nLevels = 1
triangleOptions="pAq30Dena%f" % (0.5*(L[0]*0.01*L[1]*0.01))
#level 1 : 97, 257, 161, 2.82843 
#level 2 : 354, 997, 644, 1.41421
#level 3 :1351,3926 ,2576, 0.707107 
#level 4 :5277,  15580, 10304, 0.353553 
#level 5 :20857,  62072, 41216, 0.176777 
hLevels = {0:2.82843,1:1.41421,2:0.707107,3:0.353553,4:0.176777}

#timeIntegration = BackwardEuler
#DT = 0.01*hLevels[nLevels-1]/abs(rechargeRate)#0.01
#DT = 0.01
#nDTout = int(T/DT)
#runCFL=0.1
DT = 0.001
#timeIntegrator  = testStuff.AdaptiveForwardIntegrator
#timeIntegration = testStuff.AdaptiveBackwardEuler
timeIntegration = FLCBDF
stepController  = FLCBDF_controller
systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
rtol_u[0] = 1.0e-3
atol_u[0] = 1.0e-3

timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller
systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep
runCFL = 0.1

femSpaces = {0:NC_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,2)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,1)




subgridError = None

massLumping = False

numericalFluxType = None

shockCapturing = ResGrad_SC(coefficients,nd,shockCapturingFactor=0.5,lag=True)
#shockCapturing = None

#multilevelNonlinearSolver  = NLStarILU
#multilevelNonlinearSolver  = NLGaussSeidel
#multilevelNonlinearSolver  = NLJacobi
#multilevelNonlinearSolver  = NLNI
#multilevelNonlinearSolver  = FAS
multilevelNonlinearSolver = Newton

#levelNonlinearSolver = NLStarILU
#levelNonlinearSolver = FAS
levelNonlinearSolver = Newton
#levelNonlinearSolver = NLGaussSeidel
#levelNonlinearSolver = NLJacobi

#nonlinearSmoother = NLStarILU
#nonlinearSmoother = NLGaussSeidel
nonlinearSmoother = NLJacobi

fullNewtonFlag = True

tolFac = 1.0e-6

nl_atol_res = 1.0e-6

maxNonlinearIts = 20#0#1

maxLineSearches = 10#0

matrix = SparseMatrix

multilevelLinearSolver = LU
#multilevelLinearSolver = NI

levelLinearSolver = LU
#levelLinearSolver = MGM
#numericalFluxType = Advection_DiagonalUpwind_Diffusion_IIPG_exterior

linearSmoother = Jacobi
linearSmoother = GaussSeidel
linearSmoother = StarILU

linTolFac = 0.001

conservativeFlux = {0:'p1-nc'}
