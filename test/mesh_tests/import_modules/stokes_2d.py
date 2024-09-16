from proteus import Context

##########################################################

# The default options for the context are set below.
# nLevels - this must be an iteger value >=1
# name - this must be a string
# numeric_scheme: currently implemented -
#                 TH - Taylor Hood
#                 C0P1C0P1 - stabilized continous linears velocity and pressure
# parallel - this must be a boolean value

##########################################################

opts = Context.Options([
    ("nLevels", 1, "number of levels of uniform refinement"),
    ("name", "poiseuilleFlow", "output data file name"),
    ("numeric_scheme", "TH", "specify the numerical scheme being used"),
    ("useWeakBoundaryConditions", True,"Flag: False-Strong boundary conditions, True-Weak boundary conditions"), 
    ("solveIteratively", True,"Flag: False-Direct Solver, True-Iterative Solver"),
    ("solveInParallel", False,"Flag: False-Serial Solve, True-Parallel Solve")
    ])

# TODO - add asserts

nLevels = opts.nLevels
name = opts.name
numeric_scheme = opts.numeric_scheme
useWeakBoundaryConditions = opts.useWeakBoundaryConditions
solveIteratively = opts.solveIteratively
solveInParallel = opts.solveInParallel