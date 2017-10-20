from proteus import Domain
from proteus import Context

ct=Context.Options([
    # General parameters #
    ("T",0.1,"Final time"),
    ("nDTout",1,"Number of time steps to archive"),
    ("refinement",3,"Level of refinement"),
    ("unstructured",False,"Use unstructured mesh. Set to false for periodic BCs"),
    ("SSPOrder",2,"SSP method of order 1, 2 or 3") ,
    ("cfl",0.5,"Target cfl"),
    # PARAMETERS FOR NCLS #
    ("level_set_function" ,1,"0: distance function, 1: saturated distance function"),
    ("STABILIZATION_TYPE_ncls",1,"0: SUPG, 1: EV, 2: smoothness based indicator"),
    ("SATURATED_LEVEL_SET",True,"Saturate the distance function or not?"),
    ("ENTROPY_TYPE_ncls",2,"1: parabolic, 2: logarithmic"),
    ("COUPEZ",True,"Flag to turn on/off the penalization term in Coupez approach"), 
    ("DO_REDISTANCING",True,"Solve Eikonal type equation after transport?"),    
    ("cE_ncls",1.0,"Entropy viscosity constant"),
    # PARAMETERS FOR VOF #
    ("STABILIZATION_TYPE_vof",1,"0: SUPG, 1: EV, 2: smoothness based indicator"),
    ("ENTROPY_TYPE_vof",2,"1: quadratic, 2: logarithmic"),
    ("FCT",True,"Use Flux Corrected Transport"),
    ("cE_vof",0.1,"Entropy viscosity constant"),
    ("cK",1.0,"Artificial compression constant")
],mutable=True)

# OTHER NUMERICAL PARAMETERS FOR NCLS #
epsCoupez=3
redist_tolerance=0.1
epsFactRedistance=0.33 #For signed dist function
lambda_coupez = 1.0 #Strength of redistancing and coupez force
epsFactHeaviside=epsFactDirac=1.5

# number of space dimensions #
nd=2

# MASS CORRECTION #
applyCorrection=True
correctionType = 'cg'

# General parameters #
parallel = False
linearSmoother = None
checkMass=False
runCFL = ct.cfl

# Finite elmenet spaces #
pDegree_ncls=1 
pDegree_vof=pDegree_ncls #volume of fluid should match ls for now
useHex=False
useBernstein=False

# Quadrature order #
quad_order = 2*pDegree_ncls+1

# parallel partitioning info #
from proteus import MeshTools
partitioningType = MeshTools.MeshParallelPartitioningTypes.node

# Create mesh #
nn=nnx=(2**ct.refinement)*10+1
nny=nnx
nnz=1
he=1.0/(nnx-1.0)

box=Domain.RectangularDomain(L=(1.0,1.0),
                             x=(0.0,0.0),
                             name="box");
box.writePoly("box")
if ct.unstructured:
    domain=Domain.PlanarStraightLineGraphDomain(fileprefix="box")
    domain.boundaryTags = box.boundaryTags
    bt = domain.boundaryTags
    triangleOptions="pAq30Dena%8.8f"  % (0.5*he**2,)
else:
    domain = box

# REDISTANCING #
redist_Newton=True
onlyVOF=False

# SMOOTHING FACTORS # (eps)
epsFactHeaviside=epsFactDirac=epsFact_vof=1.5
epsFactRedistance=0.33
epsFactDiffusion=10.0

# SHOCK CAPTURING PARAMETERS #
shockCapturingFactor_vof=0.2
shockCapturingFactor_ncls=0.2
shockCapturingFactor_rd=0.9
lag_shockCapturing_vof=True
lag_shockCapturing_ls=True
lag_shockCapturing_rd=False

# use absolute tolerances on al models
atolRedistance = max(1.0e-12,0.1*he)
atolConservation = max(1.0e-12,0.001*he**2)
atolVolumeOfFluid= max(1.0e-12,0.001*he**2)
atolLevelSet     = max(1.0e-12,0.001*he**2)
#controls
linearSolverConvergenceTest = 'r-true' #rits is do a set number of iterations, r-true uses true residual, PETSc default is preconditioned residual
#redist solver
fmmFlag=0

soname="cons_ls_level_"+`ct.refinement`

