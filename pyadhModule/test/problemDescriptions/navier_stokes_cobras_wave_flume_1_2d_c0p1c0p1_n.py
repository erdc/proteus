from pyadh import *
from pyadh.default_n import *
from navier_stokes_cobras_wave_flume_1_2d_p import *

#timeIntegration = FLCBDF
#stepController  = FLCBDF_controller
#systemStepControllerType = SplitOperator.Sequential_MinFLCBDFModelStep
#timeIntegration = BackwardEuler
#stepController  = FixedStep
#systemStepControllerType = SplitOperator.Sequential_FixedStep
timeIntegration = NoIntegration
stepController  = Newton_controller
rtol_u[1] = 1.0e-3
rtol_u[2] = 1.0e-3
atol_u[1] = 1.0e-2
atol_u[2] = 1.0e-2


runCFL = 10.0

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn=41
nLevels = 2
DT=1.0e-3#1.0e-3#1.0e-1#
nDTout = int(T/DT)#int(T/DT)

triangleOptions += "A"


subgridError = NavierStokesWithBodyForceASGS_velocity_pressure(coefficients,nd)


shockCapturing = None

maxNonlinearIts =2#100
maxLineSearches=10
multilevelNonlinearSolver  = NLNI#Newton#NLNI

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-4

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001

#conservativeFlux = {0:'pwl',1:'point-eval',2:'point-eval'}

