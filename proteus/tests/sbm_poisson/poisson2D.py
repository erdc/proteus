from proteus import Domain
from proteus import Context
from proteus.default_n import nLevels
from proteus.Profiling import logEvent

ct = Context.Options([
    ("degree_pol",1,"1 or 2"),
    ("parallel", False, "Use PETSc or not"),##################Use Petsc to get condition number
    ("unstructured", True, "unstructured mesh or not"),
    ("nRefine",4, "Specify initial mesh size by giving number of cells in each direction"),
    ("use_sbm",1,"0=no sbm, 1=1st order sbm, 2=2nd order sbm"),
], mutable=True)

Context.set(ct)

#===============================================================================
# Parameters
#===============================================================================
USE_SBM = ct.use_sbm

# linearSmoother = ct.linearSmoother
# compute mass balance statistics or not
# checkMass = ct.checkMass  # True
# number of space dimensions
nd = 2
# time integration, not relevant if using BDF with cfl timestepping
rtol_u = {0: 1.0e-4}
atol_u = {0: 1.0e-4}
rtol_res = {0: 1.0e-4}
atol_res = {0: 1.0e-4}
#

soname = '_'.join([str(i) for i in 
                   ['poisson_tri',`ct.degree_pol`,`ct.nRefine`,`ct.use_sbm`,`ct.unstructured`]
                   ])


parallel = ct.parallel
if parallel:
   usePETSc = True
   useSuperlu=False
else:
   usePETSc = False
#===============================================================================
# # spatial approximation orders
#===============================================================================
cDegree_ls = 0  # 0 -- CG. -1 -- DG
cDegree_vof = 0
pDegree_ls = ct.degree_pol  # level set
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
nLayersOfOverlapForParallel = 0
#===============================================================================
# spatial mesh
#===============================================================================
nn = nnx = nny = 2**ct.nRefine + 1
nnz = 1
# he = 0.6/(nnx - 1.0)
he = 0.05/2**ct.nRefine
# True for tetgen, false for tet or hex from rectangular grid
unstructured = ct.unstructured
lower_left_cornor = (-0.21, -0.21)
width_and_hight = (0.42, 0.42)
triangleFlag=1
boundaries = ['bottom', 'right', 'top', 'left', 'front', 'back', 'obstacle']
boundaryTags = dict([(key, i + 1) for (i, key) in enumerate(boundaries)])
box = Domain.RectangularDomain(L=width_and_hight,
                               x=lower_left_cornor,
                               name="box")
# box.writePoly("box")
domain = box

if unstructured:
    import my_domain
    if ct.use_sbm==1:
        domain,boundaryTags = my_domain.get_domain_2D_square(left_bottom_cornor=lower_left_cornor,
                                                       L=width_and_hight,
                                                       he=he,
                                                       )

    else:
        domain,boundaryTags = my_domain.get_domain_2D_disk(center=(0.0,0.0),
                                                       radius=0.2,
                                                       he=he,
                                                       )

    domain.writePoly("mesh")
    domain.writePLY("mesh")
#     domain.writeAsymptote("mesh")
    triangleOptions = "VApq30Dena"

#===============================================================================
# # Final time
#===============================================================================
tnList=[0,1]
#===============================================================================
# # smoothing factors
#===============================================================================
# eps
epsFactHeaviside = epsFactDirac = epsFact_vof = 1.5
epsFactRedistance = 0.33
epsFactDiffusion = 10.0
#
if useMetrics:
    shockCapturingFactor_vof = 0.5
    shockCapturingFactor_ls = 0.5
    shockCapturingFactor_rd = 0.5
    lag_shockCapturing_vof = True
    lag_shockCapturing_ls = True
    lag_shockCapturing_rd = False
else:
    shockCapturingFactor_vof = 0.2
    shockCapturingFactor_ls = 0.2
    shockCapturingFactor_rd = 0.9
    lag_shockCapturing_vof = True
    lag_shockCapturing_ls = True
    lag_shockCapturing_rd = False

# use absolute tolerances on al models
atolRedistance = max(1.0e-12, 0.1 * he)
atolConservation = max(1.0e-12, 0.001 * he**2)
atolVolumeOfFluid = max(1.0e-12, 0.001 * he**2)
atolLevelSet = max(1.0e-12, 0.001 * he**2)
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
