from proteus import Domain

#TIME INTEGRATION
timeIntegration_ls = "SSP33" #vbdf,be,flcbdf,rk
#timeIntegration_ls = "be" #vbdf,be,flcbdf,rk
runCFL = 0.3#0.3,0.185,0.125 for dgp1,dgp2,dgpk(3)
timeOrder = 1
fullNewton = False
#ENTROPY VISCOSITY
ENTROPY_VISCOSITY=1
SUPG=0
cE=5.0
cMax=0.1
#SHOCK CAPTURING
shockCapturingFactor_ls=0.0
lag_shockCapturing_ls=True

#if True uses PETSc solvers
parallel = False
linearSmoother = None
#compute mass balance statistics or not
checkMass=False#True
#number of space dimensions
nd=2
#time integration, not relevant if using BDF with cfl timestepping
rtol_u = {0:1.0e-4}
atol_u = {0:1.0e-4}
rtol_res = {0:1.0e-4}
atol_res = {0:1.0e-4}
#
#
#spatial approximation orders
cDegree_ls=0 #0 -- CG. -1 -- DG
pDegree_ls=1 #level set
useHex=False
useMetrics=0.0
vortex_quad_order = 3
#parallel partitioning info
from proteus import MeshTools
partitioningType = MeshTools.MeshParallelPartitioningTypes.node
#spatial mesh
lRefinement=3
#tag simulation name to level of refinement
#soname="vortexcgp2_bdf2_mc"+`lRefinement`
nn=nnx=nny=(2**lRefinement)*10+1
nnz=1
he=1.0/(nnx-1.0)
L=[0.0,1.0]

unstructured=False #True for tetgen, false for tet or hex from rectangular grid
box=Domain.RectangularDomain(L=(1.0,1.0),
                             x=(0.0,0.0),
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
#end time of simulation, full problem is T=8.0
T = 1.0#8.0#
#number of output time steps
nDTout = 100
#smoothing factors
#eps
epsFactHeaviside=epsFactDirac=epsFact_vof=1.5
epsFactRedistance=0.33
epsFactDiffusion=10.0
#
if useMetrics:
    shockCapturingFactor_ls=0.5
    lag_shockCapturing_ls=True

#use absolute tolerances on al models
atolRedistance = max(1.0e-12,0.1*he)
atolConservation = max(1.0e-12,0.001*he**2)
atolVolumeOfFluid= max(1.0e-12,0.001*he**2)
atolLevelSet     = max(1.0e-12,0.001*he**2)
#controls
linearSolverConvergenceTest = 'r-true' #rits is do a set number of iterations, r-true uses true residual, PETSc default is preconditioned residual
if useHex:
    hex=True
    soname="vortex_c0q"+`pDegree_ls`+"_"+timeIntegration_ls+"_level_"+`lRefinement`
else:
    soname="vortex_c0p"+`pDegree_ls`+"_"+timeIntegration_ls+"_level_"+`lRefinement`

