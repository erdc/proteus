#use special purpose solvers or not
tryOpt= True
if tryOpt:
    tryNCLS = True
    tryVOF  = True
    tryRDLS = True
    tryMCorr= True
else:
    tryNCLS = False
    tryVOF  = False
    tryRDLS = False
    tryMCorr= False
#if True uses PETSc solvers
parallel = True
#compute mass balance statistics or not
checkMass=True

#number of space dimensions
nd=2
#level set conservative or nonconservative formulation
useHJ=True
if tryNCLS:
    assert useHJ, "tryNCLS requires HJ formulation"
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
pDegree_ls=2 #level set 
cDegree_vof=0
pDegree_vof=2 #volume of fluid should match ls for now
#
#spatial quadrature orders
vortex_quad_order = 4#2*max(pDegree_vof,pDegree_ls)+1
if tryNCLS:
    if nd == 2:
        if pDegree_ls == 2:
            vortex_quad_order = 5
        else:
            vortex_quad_order = 4

    elif nd == 1:
        vortex_quad_order = 3
    else:
        if pDegree_ls == 2:
            vortex_quad_order = 5
        else:
            vortex_quad_order = 3
#parallel partitioning info
from pyadh import MeshTools
partitioningType = MeshTools.MeshParallelPartitioningTypes.node

#spatial mesh
lRefinement=3
#tag simulation name to level of refinement
soname="vortexcgp2_bdf2_mc"+`lRefinement`
if nd == 2:
    nn=nnx=nny=(2**lRefinement)*10+1;
    he=1.0/(nnx-1.0)
    L = [1.0,1.0,0.1]
elif nd == 3: 
    #cek changed to use unstructured mesh
    from tank3dDomain import *
    #this is for "2D"
    nn=nnx=nny=(2**lRefinement)*10+1;nnz=2;#11#31#26#17
    he=1.0/(nnx-1.0)
    #try to keep fixed z dimension
    L = [1.0,1.0,0.1]
    #3D should also work with just:
    #L = [1.0,1.0,1.0]
    domain = tank3d(L=L)
    bt = domain.boundaryTags
    domain.writePoly("tank3d")
    domain.writePLY("tank3d")
    domain.writeAsymptote("tank3d")
    #new flags require tetgen version >= 1.4.3
    triangleOptions="VApq1.25q12ena%21.16e" % ((he**3)/6.0,)

#end time of simulation, full problem is T=8.0
T = 8.0#8.0#
#number of output time steps
nDTout = 5
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
shockCapturingFactor_rd=0.99
#nonlinear solver absolute tolerances, all have rtol=0
#cek tightened redist tol and made conservation atol mesh-independent
#redistancing
if lRefinement >= 3:
    atolRedistance = 0.1*he
else:
    atolRedistance = 0.01*he
#mass conservation correction
atolConservation = 1.0e-6#*(he**3/6.0)
atolVolumeOfFluid= 1.0e-6
atolLevelSet     = 1.0e-6
#controls 
linearSolverConvergenceTest = 'r-true' #rits is do a set number of iterations, r-true uses true residual, PETSc default is preconditioned residual

#redist solver
fmmFlag=0
