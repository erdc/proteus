from pyadh import *
from pyadh.default_n import *
from navier_stokes_step_2d_p import *

timeIntegration = FLCBDF
stepController  = FLCBDF_controller
timeIntegration = NoIntegration
stepController  = Newton_controller
# rtol_u[1] = 1.0e-4
# rtol_u[2] = 1.0e-4
# atol_u[1] = 1.0e-4
# atol_u[2] = 1.0e-4
timeIntegration = PsiTCtte_new
tnList = [0.0,100.0*max_upstream_speed]
stepController  = PsiTCtte_controller
rtol_res[0] = 0.0
rtol_res[1] = 0.0
rtol_res[2] = 0.0
atol_res[0] = 1.0e-8
atol_res[1] = 1.0e-8
atol_res[2] = 1.0e-8



runCFL = 0.1#None
#DT=1.0e-1
DT=0.0000001
nDTout = 1#int(T/DT)

femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,2)

elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,2)

nn=3
#triangleOptions = "Aq30Dena%f" % (0.15**2 / 6.0)
triangleOptions = "Aq30Dena%f" % ((0.05*upstream_height)**2 /2.0)
triangleOptions = "Aq30Dena%f" % ((0.1*upstream_height)**2 /2.0)
#triangleOptions = "Aq30Den"
nLevels = 1

subgridError = NavierStokesASGS_velocity_pressure(coefficients,nd,lag=False)
#subgridError = StokesASGS_velocity_pressure(coefficients,nd)

massLumping = False

shockCapturing = None

numericalFluxType = None

multilevelNonlinearSolver  = NLNI
multilevelNonlinearSolver  = Newton #no nested iteration for psitc

levelNonlinearSolver = Newton
maxNonlinearIts =1000
maxLineSearches=0
tolFac = 0.0

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-5

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None
auxiliaryVariables=[RecirculationLength(upstream_length)]
