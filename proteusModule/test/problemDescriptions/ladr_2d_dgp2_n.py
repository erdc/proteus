from pyadh import *
from pyadh.default_n import *
from ladr_2d_p import *
timeOrder = 3
nStagesTime = timeOrder
runCFL = 0.01#0.1#depends on diffusion
timeIntegration = SSPRKPIintegration#LinearSSPRKintegration
#limiterType = DGlimiterP2Lagrange2d
stepController = Min_dt_RKcontroller
runCFL = 0.1#0.1
timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller

femSpaces = {0:DG_AffineQuadraticOnSimplexWithNodalBasis}
numericalFluxType = Advection_DiagonalUpwind_Diffusion_LDG#Advection_DiagonalUpwind_Diffusion_LDG#RusanovLDG


elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nnx=21; nny=21
nDTout = 101
tnList=[float(i)*T/float(nDTout-1) for i in range(nDTout)]

matrix = SparseMatrix
multilevelLinearSolver = LU

usingSSPRKNewton= False
levelNonlinearSolver = Newton#SSPRKNewton#Newton

archiveFlag = ArchiveFlags.EVERY_USER_STEP
