from __future__ import division
from past.utils import old_div
from proteus import *
from proteus.default_n import *
from proteus.mprans import SW2DCV
from . import sw_p

# ******************************************** #
# ********** READ FROM PHYSICS FILE ********** #
# ******************************************** #
nd = sw_p.nd
T = sw_p.T
nDTout = sw_p.nDTout
runCFL = sw_p.runCFL
he = sw_p.he
useSuperlu = sw_p.useSuperlu
domain = sw_p.domain
SSPOrder = sw_p.SSPOrder
LUMPED_MASS_MATRIX = sw_p.LUMPED_MASS_MATRIX
reflecting_BCs = sw_p.reflecting_BCs

# *************************************** #
# ********** MESH CONSTRUCTION ********** #
# *************************************** #
triangleFlag = sw_p.triangleFlag
nnx = sw_p.nnx
nny = sw_p.nny
nnz = 1
triangleOptions = domain.MeshOptions.triangleOptions

# ************************************** #
# ********** TIME INTEGRATION ********** #
# ************************************** #
timeIntegration = SW2DCV.RKEV
timeOrder = SSPOrder
nStagesTime = SSPOrder

# ****************************************** #
# ********** TIME STEP CONTROLLER ********** #
# ****************************************** #
stepController = Min_dt_controller

# ******************************************* #
# ********** FINITE ELEMENT SAPCES ********** #
# ******************************************* #
elementQuadrature = SimplexGaussQuadrature(nd,3)
elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis,
             1:C0_AffineLinearOnSimplexWithNodalBasis,
             2:C0_AffineLinearOnSimplexWithNodalBasis}

# ************************************** #
# ********** NONLINEAR SOLVER ********** #
# ************************************** #
multilevelNonlinearSolver  = Newton
fullNewtonFlag = False 
if (LUMPED_MASS_MATRIX==1):
    levelNonlinearSolver = ExplicitLumpedMassMatrixShallowWaterEquationsSolver
else:
    levelNonlinearSolver = ExplicitConsistentMassMatrixShallowWaterEquationsSolver
    
# ************************************ #
# ********** NUMERICAL FLUX ********** #
# ************************************ #
numericalFluxType = SW2DCV.NumericalFlux

# ************************************ #
# ********** LINEAR ALGEBRA ********** #
# ************************************ #
matrix = SparseMatrix
multilevelLinearSolver = LU
levelLinearSolver = LU
levelNonlinearSolverConvergenceTest = 'r'
linearSolverConvergenceTest = 'r-true'

# ******************************** #
# ********** TOLERANCES ********** #
# ******************************** #
nl_atol_res = 1.0e-5
nl_rtol_res = 0.0
l_atol_res = 1.0e-7
l_rtol_res = 0.0
tolFac = 0.0
maxLineSearches=0

# **************************** #
# ********** tnList ********** #
# **************************** #
tnList=[0.,1E-6]+[float(n)*T/float(nDTout) for n in range(1,nDTout+1)]
