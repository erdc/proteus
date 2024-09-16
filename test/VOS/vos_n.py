from proteus import *
from proteus.default_n import *
from vos_p import *
from vos import *
nd = 2

multilevelNonlinearSolver  = Newton
if ct.STABILIZATION_TYPE==0: #SUPG
    levelNonlinearSolver = Newton
    fullNewtonFlag = True
    updateJacobian = True
    timeIntegration = BackwardEuler_cfl
elif ct.STABILIZATION_TYPE==1:
    fullNewtonFlag = False
    updateJacobian = False
    timeIntegration = BackwardEuler_cfl
    levelNonlinearSolver = TwoStageNewton
else:
    fullNewtonFlag = False
    updateJacobian = False
    timeIntegration = VOS3P.RKEV # SSP33 
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
    if pDegree_vos == 1:
        femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
    elif pDegree_vos == 2:
        if useBernstein:
            femSpaces = {0:C0_AffineBernsteinOnCube}
        else:
            femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
    else:
        print("pDegree = %s not recognized " % pDegree_vos)
    elementQuadrature = CubeGaussQuadrature(nd,vos_quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,vos_quad_order)
else:
    if pDegree_vos == 1:
        femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
    elif pDegree_vos == 2:
        femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    else:
        print("pDegree = %s not recognized " % pDegree_vos)
    elementQuadrature = SimplexGaussQuadrature(nd,vos_quad_order)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,vos_quad_order)

#numericalFluxType = VOF.NumericalFlux
#numericalFluxType = DoNothing
numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior # PERIODIC

shockCapturing = VOS3P.ShockCapturing(coefficients,nd,shockCapturingFactor=shockCapturingFactor_vos,lag=lag_shockCapturing_vos)

matrix = SparseMatrix
if parallel:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linear_solver_options_prefix = 'vos_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU    
    levelLinearSolver = LU

if checkMass:
    auxiliaryVariables = [MassOverRegion()]

#tnList=[0.,1E-6]+[float(n)*ct.T/float(ct.nDTout) for n in range(1,ct.nDTout+1)]
tnList=[0.]+[float(n)*ct.T/float(ct.nDTout) for n in range(1,ct.nDTout+1)]

#tnList=[0.,0.0001]

#tnList=[0.,0.001,0.002,0.003]
#tnList=[0.001*i for i in range(10)]
