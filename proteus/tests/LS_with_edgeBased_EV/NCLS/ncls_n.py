from proteus import *
from proteus.default_n import *
from ncls_p import *
from ncls import *
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
    timeIntegration = NCLS.RKEV # SSP33 
    levelNonlinearSolver = ExplicitConsistentMassMatrixWithRedistancing

stepController = Min_dt_controller
runCFL = ct.cfl
timeOrder = ct.SSPOrder
nStagesTime = ct.SSPOrder

if useHex:
    hex=True
    quad=True
    if pDegree_ncls == 1:
        femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
    elif pDegree_ncls == 2:
        if useBernstein:
            femSpaces = {0:C0_AffineBernsteinOnCube}
        else:
            femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
    else:
        print "pDegree = %s not recognized " % pDegree_ncls
    elementQuadrature = CubeGaussQuadrature(nd,ncls_quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,ncls_quad_order)
else:
    if pDegree_ncls == 1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif pDegree_ncls == 2:
        femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    else:
        print "pDegree = %s not recognized " % pDegree_ncls
    elementQuadrature = SimplexGaussQuadrature(nd,ncls_quad_order)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,ncls_quad_order)

#numericalFluxType = NCLS.NumericalFlux
#numericalFluxType = DoNothing
numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior # PERIODIC

shockCapturing = NCLS.ShockCapturing(coefficients,nd,shockCapturingFactor=shockCapturingFactor_ncls,lag=lag_shockCapturing_ncls)

matrix = SparseMatrix
if parallel:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linear_solver_options_prefix = 'ncls_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU    
    levelLinearSolver = LU

if checkMass:
    auxiliaryVariables = [MassOverRegion()]

tnList=[0.,1E-6]+[float(n)*ct.T/float(ct.nDTout) for n in range(1,ct.nDTout+1)]
