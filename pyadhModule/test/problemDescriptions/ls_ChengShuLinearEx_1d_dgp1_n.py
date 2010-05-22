from pyadh import *
from pyadh.default_n import *
from ls_ChengShuLinearEx_1d_p import *


timeOrder =3
nStagesTime = timeOrder

DT=None
runCFL = 0.1#0.3#0.3#0.3 #0.1
limiterType = None
timeIntegration = SSPRKPIintegration
stepController=Min_dt_RKcontroller
usingSSPRKNewton = True

femSpaces = {0:DG_AffineLinearOnSimplexWithNodalBasis}

elementQuadrature = SimplexGaussQuadrature(nd,3)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)

nn = 41
nLevels = 3

#subgridError = None
subgridError = HamiltonJacobi_ASGS(coefficients,nd,stabFlag='1',lag=True)

massLumping = False

numericalFluxType = HamiltonJacobi_DiagonalChengShu#HamiltonJacobi_DiagonalLesaintRaviart
#HamiltonJacobi_DiagonalChengShu

shockCapturing = None
#shockCapturing = HamiltonJacobiJaffre_SC(coefficients,nd,shockCapturingFlag='1',shockCapturingFactor=0.1,lag=True)
#shockCapturing = HamiltonJacobi_SC(coefficients,nd,shockCapturingFlag='2',shockCapturingFactor=0.05,lag=True)#None



multilevelNonlinearSolver  = NLNI
maxNonlinearIts = 20
levelNonlinearSolver = SSPRKNewton#Newton

nonlinearSmoother = NLGaussSeidel

fullNewtonFlag = True

tolFac = 0.01

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU

levelLinearSolver = LU

linearSmoother = GaussSeidel

linTolFac = 0.001

conservativeFlux = None


archiveFlag = ArchiveFlags.EVERY_USER_STEP
