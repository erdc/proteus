#if True uses PETSc solvers
import os 

parallel = False
linearSmoother = None
#compute mass balance statistics or not
checkMass=False#True
#number of space dimensions
nd=3
#time integration, not relevant if using BDF with cfl timestepping
rtol_u = {0:1.0e-4}
atol_u = {0:1.0e-4}
rtol_res = {0:1.0e-4}
atol_res = {0:1.0e-4}
#
timeIntegration_vof = "BE"
timeIntegration_ls = "BE"
#if want bdf2 or bdf1
timeOrder = 2
runCFL = 0.3#0.3,0.185,0.125 for dgp1,dgp2,dgpk(3)
#
#spatial approximation orders
cDegree_ls=0 #0 -- CG. -1 -- DG
cDegree_vof=0
pDegree_ls=1 #level set
pDegree_vof=pDegree_ls #volume of fluid should match ls for now
useHex=False#True
#
#spatial quadrature orders
#2*max(pDegree_vof,pDegree_ls)+1
if pDegree_ls == 2:
    vortex_quad_order = 5
else:
    vortex_quad_order = 3
#parallel partitioning info
from proteus import MeshTools
partitioningType = MeshTools.MeshParallelPartitioningTypes.node
#spatial mesh
lRefinement=1
#tag simulation name to level of refinement
#soname="vortexcgp2_bdf2_mc"+`lRefinement`
pseudo2D=True
if pseudo2D:
    nn=nnx=nny=(2**lRefinement)*5+1
    nnz=2
    he = 0.5
    #he=1.0/(nnx-1.0)
    L=[1.0,1.0,he]
else:
    nn=nnx=nny=nnz=(2**lRefinement)*10+1
    he = 1.0/(nnx-1.0)
    L = [1.0,1.0,1.0]
unstructured=True#True for tetgen, false for tet or hex from rectangular grid
genMesh=False
if unstructured:
    from .tank3dDomain import *
    domain = tank3d(L=L)
    bt = domain.boundaryTags
    #domain.writePoly("tank3d")
    #domain.writePLY("tank3d")
    #domain.writeAsymptote("tank3d")
    domain.polyfile=os.path.dirname(os.path.abspath(__file__))+"/"+"tank3d"
    triangleOptions="VApq1.3q18ena%21.16e" % ((he**3)/6.0,)
else:
    from proteus.Domain import RectangularDomain
    domain = RectangularDomain(L)
#end time of simulation, full problem is T=8.0
T = 8.0#8.0#
#number of output time steps
nDTout = 80
#mass correction
applyCorrection=True
applyRedistancing=True
#smoothing factors
#eps
epsFactHeaviside=1.5
epsFactDirac=1.5
epsFactDiffusion=10.0
epsFactRedistance=0.33
epsFact_vof=1.5
#
shockCapturingFactor_vof=0.33
shockCapturingFactor_ls=0.33
shockCapturingFactor_rd=0.9
#use absolute tolerances on al models
atolRedistance = 1.0e-3
atolConservation = 1.0e-6
atolVolumeOfFluid= 1.0e-4
atolLevelSet     = 1.0e-4
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
    soname="vortex_c0q"+repr(pDegree_ls)+correctionType+"_bdf_"+repr(timeOrder)+"_level_"+repr(lRefinement)
else:
    soname="vortex_c0p"+repr(pDegree_ls)+correctionType+"_bdf_"+repr(timeOrder)+"_level_"+repr(lRefinement)
