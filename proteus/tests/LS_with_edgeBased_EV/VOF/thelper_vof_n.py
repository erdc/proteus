from __future__ import print_function
from __future__ import absolute_import
from builtins import range
from proteus import *
from proteus.default_n import *
from .thelper_vof_p import *
from .thelper_vof import *
nd = 2

multilevelNonlinearSolver  = Newton
if ct.STABILIZATION_TYPE==0: #SUPG
    levelNonlinearSolver = Newton
    fullNewtonFlag = True
    updateJacobian = True
    timeIntegration = BackwardEuler_cfl
else:
    fullNewtonFlag = False
    updateJacobian = False
    timeIntegration = VOF.RKEV # SSP33 
    if ct.LUMPED_MASS_MATRIX==True: 
        levelNonlinearSolver = ExplicitLumpedMassMatrix
    else:
        levelNonlinearSolver = ExplicitConsistentMassMatrixForVOF

stepController = Min_dt_controller
runCFL = ct.cfl
timeOrder = ct.SSPOrder
nStagesTime = ct.SSPOrder

if useHex:
    hex=True
    quad=True
    if pDegree_vof == 1:
        femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
    elif pDegree_vof == 2:
        if useBernstein:
            femSpaces = {0:C0_AffineBernsteinOnCube}
        else:
            femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
    else:
        print("pDegree = %s not recognized " % pDegree_vof)
    elementQuadrature = CubeGaussQuadrature(nd,vof_quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,vof_quad_order)
else:
    if pDegree_vof == 1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif pDegree_vof == 2:
        femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    else:
        print("pDegree = %s not recognized " % pDegree_vof)
    elementQuadrature = SimplexGaussQuadrature(nd,vof_quad_order)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vof_quad_order)

#numericalFluxType = VOF.NumericalFlux
#numericalFluxType = DoNothing
numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior # PERIODIC

shockCapturing = VOF.ShockCapturing(coefficients,nd,shockCapturingFactor=shockCapturingFactor_vof,lag=lag_shockCapturing_vof)

matrix = SparseMatrix
if parallel:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linear_solver_options_prefix = 'vof_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU    
    levelLinearSolver = LU

if checkMass:
    auxiliaryVariables = [MassOverRegion()]

tnList=[0.,1E-6]+[float(n)*ct.T/float(ct.nDTout) for n in range(1,ct.nDTout+1)]

