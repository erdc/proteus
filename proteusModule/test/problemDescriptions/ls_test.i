# uncomment the next  two lines to use dense matrices instead  of sparse
#import Numeric
#n.matrix = Numeric.array
n.multilevelLinearSolver = LinearSolvers.LU
n.levelLinearSolver = LinearSolvers.LU
start
n.multilevelLinearSolver = LinearSolvers.GaussSeidel
n.levelLinearSolver = LinearSolvers.GaussSeidel
start
n.multilevelLinearSolver = LinearSolvers.NI
n.levelLinearSolver = LinearSolvers.GaussSeidel
start
n.multilevelLinearSolver = LinearSolvers.NI
n.levelLinearSolver = LinearSolvers.MGM
n.linearSmoother = LinearSolvers.GaussSeidel
start
# uncomment the next lines to test jacobi, nested jacobi, and multigrid with jacobi smoothing
# n.multilevelLinearSolver = LinearSolvers.Jacobi
# n.levelLinearSolver = LinearSolvers.Jacobi
# s
# n.multilevelLinearSolver = LinearSolvers.NI
# n.levelLinearSolver = LinearSolvers.Jacobi
# s
# n.multilevelLinearSolver = LinearSolvers.NI
# n.levelLinearSolver = LinearSolvers.MGM
# n.linearSmoother = LinearSolvers.Jacobi
# s
# uncomment the next lines to test ASM
# n.multilevelLinearSolver = LinearSolvers.StarILU
# n.levelLinearSolver = LinearSolvers.StarILU
# s
# n.multilevelLinearSolver = LinearSolvers.NI
# n.levelLinearSolver = LinearSolvers.StarILU
# s
# n.multilevelLinearSolver = LinearSolvers.NI
# n.levelLinearSolver = LinearSolvers.MGM
# n.linearSmoother = LinearSolvers.StarILU
# s
quit

