from proteus import Domain
from proteus import Context
from proteus.default_n import nLevels

ct = Context.Options([
    ("T", 8.0, "Time interval [0, T]"),
    ("checkMass", False, "Check mass or not"),
    ("cfl", 0.3, "Target CFL number"),
    ("parallel", False, "Use PETSc or not"),
    ("lRefinement", 3, "Specify initial mesh size by giving number of cells in each direction"),
    ("applyCorrection", True, "apply mass correction or not"),
    ("applyRedistancing", True, "apply redkstancing or not"),
    ("hk", 1.0, "the size of sub-triangles in the reference triangle used for comp quadrature rule"),
    ("massCorrectionReference", 4, "reference mass in MCorr"),
    ("thetaMCorr", 0.5, " >=0, <=1, e.g., 0.0=explicit, 0.5=Crank-Nicolson, 1.0=implicit"),
    ("degreeFE",1,"degree of Lagrangian polynomial on each element")
], mutable=True)


#if True uses PETSc solvers
parallel = ct.parallel
linearSmoother = None
#compute mass balance statistics or not
checkMass=ct.checkMass#True
#number of space dimensions
nd=2
#time integration, not relevant if using BDF with cfl timestepping
rtol_u = {0:1.0e-4}
atol_u = {0:1.0e-4}
rtol_res = {0:1.0e-4}
atol_res = {0:1.0e-4}
#
timeIntegration_vof = "vbdf"#vbdf,be,flcbdf,rk
timeIntegration_ls = "vbdf"#vbdf,be,flcbdf,rk
timeOrder = 2

runCFL = ct.cfl#0.3,0.185,0.125 for dgp1,dgp2,dgpk(3)
#
#spatial approximation orders
cDegree_ls=0 #0 -- CG. -1 -- DG
cDegree_vof=0
pDegree_ls=ct.degreeFE #level set
pDegree_vof=pDegree_ls #volume of fluid should match ls for now
useHex=False#True
useMetrics=1.0
#
#spatial quadrature orders
#2*max(pDegree_vof,pDegree_ls)+1
if pDegree_ls == 2:
    vortex_quad_order = 4
else:
    vortex_quad_order = 3

# mass correction
mass_correction_reference = ct.massCorrectionReference #: 0=H(phi^{n+1}); 1=explicit-1; 2=implicit; 3=explicit-2; 4=theta method
theta_time_discretization_mcorr = ct.thetaMCorr#: if mass_correction_reference==4

#sub-element edge size, used to create composite quadrature rule
hk = ct.hk

#parallel partitioning info
from proteus import MeshTools
partitioningType = MeshTools.MeshParallelPartitioningTypes.node
#spatial mesh
lRefinement=ct.lRefinement
#tag simulation name to level of refinement
#soname="vortexcgp2_bdf2_mc"+`lRefinement`
nn=nnx=nny=(2**lRefinement)*10+1
nnz=1
he=1.0/(nnx-1.0)
L=[1.0,1.0]

unstructured=False#True for tetgen, false for tet or hex from rectangular grid
if unstructured:
    from tank2dDomain import *
    domain = tank2d(L=L)
    bt = domain.boundaryTags
    domain.writePoly("tank2d")
    triangleOptions="pAq30Dena%8.8f"  % (0.5*he**2,)
else:
    from proteus.Domain import RectangularDomain
    domain = RectangularDomain(L)
#end time of simulation, full problem is T=8.0
T = ct.T#8.0#
#number of output time steps
nDTout = 80
#mass correction
applyCorrection=ct.applyCorrection
applyRedistancing=ct.applyRedistancing
redist_Newton=True
onlyVOF=False#True
#smoothing factors
#eps
epsFactHeaviside=epsFactDirac=epsFact_vof=1.5*hk
epsFactRedistance=0.33
epsFactDiffusion=10.0
#
if useMetrics:
    shockCapturingFactor_vof=0.25
    shockCapturingFactor_ls=0.25
    shockCapturingFactor_rd=0.5
    lag_shockCapturing_vof=True
    lag_shockCapturing_ls=True
    lag_shockCapturing_rd=False
else:
    shockCapturingFactor_vof=0.2
    shockCapturingFactor_ls=0.2
    shockCapturingFactor_rd=0.9
    lag_shockCapturing_vof=True
    lag_shockCapturing_ls=True
    lag_shockCapturing_rd=False

#use absolute tolerances on al models
atolRedistance = max(1.0e-12,0.01*he)
#atolConservation = max(1.0e-12,0.001*he**2)#1.0e-12#
#atolVolumeOfFluid= max(1.0e-12,0.001*he**2)
#atolLevelSet     = max(1.0e-12,0.001*he**2)
atolConservation = 1.0e-10
atolVolumeOfFluid= 1.0e-10
atolLevelSet     = 1.0e-10
#controls
linearSolverConvergenceTest = 'r-true' #rits is do a set number of iterations, r-true uses true residual, PETSc default is preconditioned residual
#redist solver
fmmFlag=0
#
#correctionType = 'dg'
#correctionType = 'dgp0'
#correctionType = 'global'
correctionType = 'cg'
#correctionType = 'none'
if useHex:
    hex=True
    soname="vortex_c0q"+`pDegree_ls`+correctionType+"_"+timeIntegration_vof+"_"+`timeOrder`+"_level_"+`lRefinement`
else:
    soname="vortex_c0p"+`pDegree_ls`+correctionType+"_"+timeIntegration_vof+"_"+`timeOrder`+"_level_"+`lRefinement`+"_"+`ct.applyCorrection`+"_"+`ct.applyRedistancing`+"_"+`hk`+"_"+`mass_correction_reference`+"_"+`ct.thetaMCorr`
