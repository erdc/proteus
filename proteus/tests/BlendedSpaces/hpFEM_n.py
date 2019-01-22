from proteus import *
from proteus.default_n import *
from hpFEM_p import *
from hpFEM import *
nd = 2

##########
# SOLVER #
##########
multilevelNonlinearSolver  = Newton
levelNonlinearSolver = myLinearSolver
fullNewtonFlag = True
updateJacobian = True

#######################
# ABOUT TIME STEPPING # This is dummy
#######################
timeIntegration = BackwardEuler_cfl
stepController = Min_dt_controller
runCFL = 0.33 # this is dummy

###################
# QUADRATURE RULE #
###################
if COMPOSITE_QUAD_RULE:
    quad_order = 4
else:
    quad_order = 2*(pDegree)+2
#

#########################
# FINITE ELEMENT SPACES #
#########################
if useHex==True:
    hex=True
    quad=True
    if pDegree == 1:
        femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
        auxSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
    elif pDegree == 2:
        auxSpaces = {0:C0_AffineLinearOnQuadraticCube}
        if useLinearOnQuadraticCube:
            femSpaces = {0:C0_AffineLinearOnQuadraticCube}
        elif useBernstein:
            femSpaces = {0:C0_AffineBernsteinOnCube}
        else:
            femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
    else:
        print "pDegree = %s not recognized " % pDegree
    #
    if COMPOSITE_QUAD_RULE:
        elementQuadrature = CubeGaussCompositeQuadrature(nd,quad_order)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,quad_order*2)
    else:
        elementQuadrature = CubeGaussQuadrature(nd,quad_order)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,quad_order)
else:
    raise NotImplementedError("Only works with quads")

numericalFluxType = DoNothing

#####################
# FOR PARALLEL RUNS # NOTE: the code is not tested in parallel
#####################
matrix = SparseMatrix
if parallel:
    multilevelLinearSolver = KSP_petsc4py
    levelLinearSolver = KSP_petsc4py
    linear_solver_options_prefix = 'poisson_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

tnList=[0.,0.1] #this is dummy
