from pyadh.LinearSolvers  import *
#multigrid on
n.multilevelLinearSolver = NI
n.levelLinearSolver = MGM
n.smoother = StarILU
start
n.smoother = GaussSeidel
start
n.smoother = Jacobi
start
#multigrid off
n.multilevelLinearSolver = NI
n.levelLinearSolver = StarILU
start
n.levelLinearSolver = GaussSeidel
start
n.levelLinearSolver = Jacobi
start
#nested iteration off
n.multilevelLinearSolver = StarILU
start
n.multilevelLinearSolver = GaussSeidel
start
n.multilevelLinearSolver = Jacobi
start
quit

