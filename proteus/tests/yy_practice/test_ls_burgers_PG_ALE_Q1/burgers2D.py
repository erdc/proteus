from proteus import Domain
from proteus import Context
from proteus.default_n import nLevels

ct = Context.Options([
    ("T", 0.5, "Time interval [0, T]"),
    ("nDTout", 40, "Times to print out [0,T]"),
    ("checkMass", False, "Check mass or not"),
    ("cfl", 0.3, "Target CFL number"),
    ("parallel", False, "Use PETSc or not"),
    ("linearSmoother", False, "Use linear smoother or not"),
    ("unstructured", False, "unstructured mesh or not"),
    ("ncells", 32, "Specify initial mesh size by giving number of cells in each direction"),
    ("nLevels", 1, "number of refiments"),
    ("timeIntegration_ls", 'be',
     "method for Time integration: 'be', 'vbdf', 'flcbdf','rk' "),
    ("datafile", 'errorInfo.db', "Filename to save error information"),
    ("useHex", True, "use quadrilateral or not"),
    ("stablization", 1, "Stabilization method: 0=SUPG, 1=EV, 2=FCT"),
    ("method", 1, "Method used to solve ALE: 0=midified ML, 1=piecewise linear moving vertex")
], mutable=True)

if ct.useHex:
    hex = True
    quad = True
    soname = '_'.join(
        [str(i) for i in ['burgers_qua', ct.timeIntegration_ls, ct.ncells, ct.nLevels]])
else:
    soname = '_'.join(
        [str(i) for i in ['burgers_tri', ct.timeIntegration_ls, ct.ncells, ct.nLevels]])

ct.datafile = soname + '_' + ct.datafile

Context.set(ct)


# if True uses PETSc solvers
parallel = ct.parallel
linearSmoother = ct.linearSmoother
# compute mass balance statistics or not
checkMass = ct.checkMass  # True
# number of space dimensions
nd = 2
# time integration, not relevant if using BDF with cfl timestepping
rtol_u = {0: 1.0e-4}
atol_u = {0: 1.0e-4}
rtol_res = {0: 1.0e-4}
atol_res = {0: 1.0e-4}
#


runCFL = ct.cfl  # 0.3,0.185,0.125 for dgp1,dgp2,dgpk(3)
#
# spatial approximation orders
cDegree_ls = 0  # 0 -- CG. -1 -- DG
cDegree_vof = 0
pDegree_ls = 1  # level set
pDegree_vof = pDegree_ls  # volume of fluid should match ls for now
useHex = False  # True
useMetrics = 1.0
#
# spatial quadrature orders
# 2*max(pDegree_vof,pDegree_ls)+1
if pDegree_ls == 2:
    rotation_quad_order = 5
else:
    rotation_quad_order = 3
# parallel partitioning info
from proteus import MeshTools
partitioningType = MeshTools.MeshParallelPartitioningTypes.node
# spatial mesh
# lRefinement = 1
# tag simulation name to level of refinement
# soname="rotationcgp2_bdf2_mc"+`lRefinement`
nn = nnx = nny = ct.ncells + 1
nnz = 1
he = 1.0 / (nnx - 1.0)

# True for tetgen, false for tet or hex from rectangular grid
unstructured = ct.unstructured

L = lower_left_cornor = (-0.2, -.2)
width_and_hight = (2.0, 2.0)
rotation_center = (lower_left_cornor[0] + 0.5 * width_and_hight[0],
                   lower_left_cornor[1] + 0.5 * width_and_hight[1])

initial_area = width_and_hight[0] * width_and_hight[1] / nnx / nny / 2.0

box = Domain.RectangularDomain(L=width_and_hight,
                               x=lower_left_cornor,
                               name="box")
box.writePoly("box")
if unstructured:
    from rotationDomain import *
    domain = Domain.PlanarStraightLineGraphDomain(fileprefix="box")
    domain.boundaryTags = box.boundaryTags
    bt = domain.boundaryTags
    triangleOptions = "pAq30Dena%8.8f" % (0.5 * he**2,)
else:
    domain = box

# Final time
T = ct.T
# number of output time steps
nDTout = ct.nDTout
# mass correction
# applyCorrection = False
# applyRedistancing = False
# redist_Newton = False
# onlyVOF = False  # True
# smoothing factors
# eps
epsFactHeaviside = epsFactDirac = epsFact_vof = 1.5
epsFactRedistance = 0.33
epsFactDiffusion = 10.0
#
# if useMetrics:
#     shockCapturingFactor_vof = 0.5
#     shockCapturingFactor_ls = 0.5
#     shockCapturingFactor_rd = 0.5
#     lag_shockCapturing_vof = True
#     lag_shockCapturing_ls = True
#     lag_shockCapturing_rd = False
# else:
#     shockCapturingFactor_vof = 0.2
#     shockCapturingFactor_ls = 0.2
#     shockCapturingFactor_rd = 0.9
#     lag_shockCapturing_vof = True
#     lag_shockCapturing_ls = True
#     lag_shockCapturing_rd = False

# # use absolute tolerances on al models
# atolRedistance = max(1.0e-12, 0.1 * he)
# atolConservation = max(1.0e-12, 0.001 * he**2)
# atolVolumeOfFluid = max(1.0e-12, 0.001 * he**2)
# atolLevelSet = max(1.0e-12, 0.001 * he**2)
# controls
# rits is do a set number of iterations, r-true uses true residual, PETSc
# default is preconditioned residual
linearSolverConvergenceTest = 'r-true'
# redist solver
fmmFlag = 0
#
#correctionType = 'dg'
#correctionType = 'dgp0'
#correctionType = 'global'
# correctionType = ct.correctionType
#correctionType = 'none'
# if ct.useHex:
#     hex = True
#     soname = "rotation_qua_c0q" + `pDegree_ls`+"_" + \
#         ct.timeIntegration_vof + "_level_" + `ct.nLevels`
# else:
#     soname = "rotation_tri_c0p" + `pDegree_ls`+"_" + \
#         ct.timeIntegration_vof + "_level_" + `ct.nLevels`
