from proteus import Context

# TODO - add weak/strong and direct/indirect solver options

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
    ("name", "drivenCavityNSETrial", "output data file name"),
    ("numeric_scheme", "THQuads", "specify the numerical scheme being used"),
    ("useWeakBoundaryConditions", True, "Flag: False-Strong boundary conditions, True-Weak boundary conditions"),
    ("solveIteratively", True, "Flag: False-Direct Solver, True-Iterative Solver"),
    ("solveInParallel", False,"Flag: False-Serial Solve, True-Parallel Solve"),
    ("schur_solver", "Qp", "Flag: Options - selfp, Qp, PCD, LSC"),
    ("RE_through_bdy",True,"Flag: is the reynolds number enforced through boundary condition?")
    ])

# TODO - add asserts and documentation explaining options

nLevels = opts.nLevels
name = opts.name
numeric_scheme = opts.numeric_scheme
useWeakBoundaryConditions = opts.useWeakBoundaryConditions
solveIteratively = opts.solveIteratively
solveInParallel = opts.solveInParallel
schur_solver = opts.schur_solver
RE_through_bdy = opts.RE_through_bdy