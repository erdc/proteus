from pyadh import *
from pyadh.default_n import *
from sw_dam_break_santillana_2d_p import *

nLevels=1
triangleOptions="q30Dena%f" % (0.5*(0.1**2),)
tnList=[0.0,dt_init,T]
nDTout = 101
tnList = [0.0 + i*T/(float(nDTout)-1.) for i in range(nDTout)]
archiveFlag = ArchiveFlags.EVERY_USER_STEP

runCFL=0.45#0.15
timeOrder = 1
timeIntegration = LinearSSPRKintegration 
stepController=Min_dt_RKcontroller
systemStepControllerType = SplitOperator.Sequential_MinAdaptiveModelStep

femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis,
             1:DG_AffineP0_OnSimplexWithMonomialBasis,
             2:DG_AffineP0_OnSimplexWithMonomialBasis}

elementQuadrature = SimplexGaussQuadrature(nd,1)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,1)
numericalFluxType = ShallowWater_2D
multilevelNonlinearSolver  = NLNI
usingSSPRKNewton=True
levelNonlinearSolver = SSPRKNewton#Newton
fullNewtonFlag = True

subgridError = ShallowWater_CFL(coefficients,nd,g)

massLumping=False

shockCapturing = None



tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix
#matrix = numpy.array
multilevelLinearSolver = LU

levelLinearSolver = LU

linTolFac = 0.001

#conservativeFlux = {0:'pwl'}
parallelPartitioningType = MeshParallelPartitioningTypes.node

archiveFlag = ArchiveFlags.EVERY_USER_STEP
