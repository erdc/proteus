from proteus import *
from proteus.default_n import *
from clsvof_p import *
from clsvof import *
nd = 3

# ABOUT NONLINEAR SOLVER #
multilevelNonlinearSolver  = Newton
#levelNonlinearSolver = Newton
levelNonlinearSolver = TwoStageNewton_for_CLSVOF
fullNewtonFlag = True
updateJacobian = True
nl_atol_res=1E-12
tolFac=0.
maxNonlinearIts = 100000

# ABOUT TIME INTEGRATION (BASE) #
timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller
runCFL = cfl

# FINITE ELEMENT SPACE #
if useHex:
    hex=True
    quad=True
    if pDegree_clsvof == 1:
        femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
    elif pDegree_clsvof == 2:
        if useBernstein:
            femSpaces = {0:C0_AffineBernsteinOnCube}
        else:
            femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
    else:
        print "pDegree = %s not recognized " % pDegree_clsvof
    elementQuadrature = CubeGaussQuadrature(nd,clsvof_quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,clsvof_quad_order)
else:
    if pDegree_clsvof == 1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif pDegree_clsvof == 2:
        if useBernstein:
            femSpaces = {0:C0_AffineBernsteinOnSimplex}
        else:                
            femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    else:
        print "pDegree = %s not recognized " % pDegree_clsvof
    elementQuadrature = SimplexGaussQuadrature(nd,clsvof_quad_order)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,clsvof_quad_order)

#numericalFluxType = CLSVOF.NumericalFlux
#numericalFluxType = DoNothing
numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior # PERIODIC

matrix = SparseMatrix
if parallel:
    multilevelLinearSolver = PETSc
    levelLinearSolver = PETSc
    #multilevelLinearSolver = KSP_petsc4py
    #levelLinearSolver = KSP_petsc4py
    linear_solver_options_prefix = 'clsvof_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU    
    levelLinearSolver = LU

#tnList=[0.,1E-6]+[float(n)*ct.T/float(ct.nDTout) for n in range(1,ct.nDTout+1)]
tnList=[float(n)*ct.T/float(ct.nDTout) for n in range(0,ct.nDTout+1)]

