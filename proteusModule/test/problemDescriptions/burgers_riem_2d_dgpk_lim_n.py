from pyadh import *
from pyadh.default_n import *
from burgers_riem_2d_p import *

timeOrder =3
nStagesTime = timeOrder

DT=None
runCFL = 0.185 #0.1

limiterType = TimeIntegration.DGlimiterPkMonomial2d#

timeIntegration = SSPRKPIintegration
stepController=Min_dt_RKcontroller
nDTout = 1

femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}

elementQuadrature = SimplexGaussQuadrature(nd,4)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)

nn = 11
nLevels = 1

subgridError = None

massLumping = False

class BurgersNumericalFlux(ConvexOneSonicPointNumericalFlux):
   def __init__(self,vt,getPointwiseBoundaryConditions,
                 getAdvectiveFluxBoundaryConditions,
                 getDiffusiveFluxBoundaryConditions):
        ConvexOneSonicPointNumericalFlux.__init__(self,vt,getPointwiseBoundaryConditions,
                                                  getAdvectiveFluxBoundaryConditions,
                                                  getDiffusiveFluxBoundaryConditions,
                                                  sonicPoint=0.0,sonicFlux=0.0)
   #
#
numericalFluxType =  BurgersNumericalFlux#RusanovNumericalFlux_Diagonal#Advection_DiagonalUpwind#BurgersNumericalFlux#

shockCapturing = None

multilevelNonlinearSolver  = NLNI
usingSSPRKNewton=True
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

checkMass = True

archiveFlag = ArchiveFlags.EVERY_USER_STEP
