from proteus import Domain

timeIntegration_vof = "SSP33" #vbdf,be,flcbdf,rk
#timeIntegration_vof = "be" #vbdf,be,flcbdf,rk
fullNewton=False
#ENTROPY VISCOSITY and ART COMPRESSION PARAMETERS
ENTROPY_VISCOSITY=1
SUPG=0
cE = 10.0
cMax = 0.1
cK = 1.0
shockCapturingFactor_vof=0.2
#Other time parameters
timeOrder = 1
runCFL = 0.1#0.3,0.185,0.125 for dgp1,dgp2,dgpk(3)
lag_shockCapturing_vof=True
#if True uses PETSc solvers
parallel = False
linearSmoother = None
#compute mass balance statistics or not
checkMass=False
#number of space dimensions
nd=2
#time integration, not relevant if using BDF with cfl timestepping
rtol_u = {0:1.0e-4}
atol_u = {0:1.0e-4}
rtol_res = {0:1.0e-4}
atol_res = {0:1.0e-4}
#
#spatial approximation orders
cDegree_vof=0
pDegree_vof=1
useHex=False
useMetrics=0.0
#
#spatial quadrature orders
rotation_quad_order = 2*pDegree_vof+1
#parallel partitioning info
from proteus import MeshTools
partitioningType = MeshTools.MeshParallelPartitioningTypes.node
#spatial mesh
lRefinement=3
#tag simulation name to level of refinement
nn=nnx=nny=(2**lRefinement)*10+1
nnz=1
he=1.0/(nnx-1.0)
L=[2.0,2.0]

unstructured=False #True for tetgen, false for tet or hex from rectangular grid
box=Domain.RectangularDomain(L=(2.0,2.0),
                             x=(-1.0,-1.0),
                             name="box");
box.writePoly("box")
if unstructured:
    from tank2dDomain import *
    domain=Domain.PlanarStraightLineGraphDomain(fileprefix="box")
    domain.boundaryTags = box.boundaryTags
    bt = domain.boundaryTags
    triangleOptions="pAq30Dena%8.8f"  % (0.5*he**2,)
else:
    domain = box
#end time of simulation
T = 1.0
#number of output time steps
nDTout = 10
#smoothing factors
#eps
epsFactHeaviside=epsFactDirac=epsFact_vof=1.5 #1.5
#
if useMetrics:
    shockCapturingFactor_vof=0.5
    lag_shockCapturing_vof=True

#use absolute tolerances on al models
atolRedistance = max(1.0e-12,0.1*he)
atolConservation = max(1.0e-12,0.001*he**2)
atolVolumeOfFluid= max(1.0e-12,0.001*he**2)
atolLevelSet     = max(1.0e-12,0.001*he**2)
#controls
linearSolverConvergenceTest = 'r-true' #rits is do a set number of iterations, r-true uses true residual, PETSc default is preconditioned residual
#redist solver
fmmFlag=0
#
if useHex:
    hex=True
    soname="rotation_c0q"+`pDegree_vof`+"_"+timeIntegration_vof+"_level_"+`lRefinement`
else:
    soname="rotation_c0p"+`pDegree_vof`+"_"+timeIntegration_vof+"_level_"+`lRefinement`
