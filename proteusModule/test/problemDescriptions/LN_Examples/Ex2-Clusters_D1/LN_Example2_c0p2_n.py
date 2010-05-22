from pyadh import *
from pyadh.default_n import *
from LN_Example2_p import *

timeIntegration = NoIntegration
nDTout = 1

femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}

subgridError = None

shockCapturing = None

multilevelNonlinearSolver  = Newton

levelNonlinearSolver = Newton

fullNewtonFlag = True

tolFac = 0.0

nl_atol_res = 1.0e-8

matrix = SparseMatrix

multilevelLinearSolver = LU
#multilevelLinearSolver = PETSc

levelLinearSolver = LU
#levelLinearSolver = PETSc


l_atol_res = 1.0e-10
linTolFac = 0.0

numericalFluxType = numericalFlux_flow

conservativeFlux = {0:vpp_flag}#{0:'pwl-ib-fix-0'}#{0:'pwl'}{0:'pwc'}

#calculate velocity norms over region with low permeability
auxiliaryVariables = [AuxiliaryVariables.VelocityNormOverRegion(regionIdList=[1,2,3,4,5,6])]
