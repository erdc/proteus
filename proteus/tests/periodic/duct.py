import numpy as np
from proteus.defaults import (Physics_base,
                              Numerics_base,
                              System_base)
from proteus import (Domain,
                     TimeIntegration,
                     StepControl,
                     FemTools,
                     Quadrature,
                     NonlinearSolvers,
                     LinearAlgebraTools,
                     LinearSolvers,
                     MeshTools,
                     Context,
                     AnalyticalSolutions)
from proteus.default_so import *
from proteus.mprans import RANS2P

opts = Context.Options([
    ("nd", 2, "Number of space dimensions"),
    ("grid", True, "Use a regular grid"),
    ("triangles", True, "Use triangular or tetrahedral elements"),
    ("spaceOrder", 1, "Use (bi-)linear or (bi-)quadratic spaces"),
    ("timeOrder", 1, "Use bdf1 or bdf2"),
    ("periodic", False, "Use periodic boundary conditions"),
    ("weak", True, "Use weak boundary conditions"),
    ("coord", False, "Use coordinates for setting boundary conditions"),
    ("Re", 1.0, "Reynolds number for non-periodic problem"),
    ("nnx", 21, "Number of grid nodes in x-direction"),
    ("pc_type", 'LU', "Specify preconditioner type"),
    ("A_block_AMG", False, "Specify whether a block-AMG should be used to solve A-block")])

#########
#Physics#
#########
p = Physics_base(nd = opts.nd,
                 name="duct_physics")

p.L=(4.0, 1.0, 1.0)
if p.nd == 2:
    p.L=(4.0, 1.0)
p.domain = Domain.RectangularDomain(p.L)
boundaryTags = p.domain.boundaryTags
if not opts.grid:
    p.domain.writePoly("duct")
    if p.nd == 3:
        p.domain = Domain.PiecewiseLinearComplexDomain()
    elif p.nd == 2:
        p.domain = Domain.PlanarStraightLineGraphDomain()
    p.domain.readPoly("duct")

nu = 1.004e-6
rho = 998.2

if p.nd == 3:
    h = p.L[2]
    umax = opts.Re*nu/p.L[2]
else:
    h = p.L[1]
    umax = opts.Re*nu/p.L[1]
p.T = 2.0*(p.L[0]/umax)
mu = nu*rho
G = (umax*8.0*mu)/(h**2)

if opts.periodic:
    gravity = [G/rho, 0., 0.]
else:
    gravity = [0., 0., 0.]

p.LevelModelType = RANS2P.LevelModel

if opts.periodic:
    nullSpace="NavierStokesConstantPressure"
else:
    nullSpace="NoNullSpace"


p.coefficients = RANS2P.Coefficients(epsFact=0.0,
                                     sigma=0.0,
                                     rho_0=rho,nu_0=nu,
                                     rho_1=rho,nu_1=nu,
                                     g=gravity,
                                     nd=p.nd,
                                     ME_model=0,
                                     VF_model=None,
                                     LS_model=None,
                                     Closure_0_model=None,
                                     Closure_1_model=None,
                                     epsFact_density=0.0,
                                     stokes=False,
                                     useVF=0.0,
                                     useRBLES=0.0,
                                     useMetrics=1.0,
                                     eb_adjoint_sigma=1.0,
                                     eb_penalty_constant=100.0,
                                     forceStrongDirichlet=not opts.weak,
                                     turbulenceClosureModel=0,
                                     NONCONSERVATIVE_FORM=1.0,
                                     MOMENTUM_SGE=1.0,
                                     PRESSURE_SGE=1.0,
                                     VELOCITY_SGE=1.0,
                                     nullSpace=nullSpace)

eps=1.0e-8
if opts.periodic:
    if opts.nd == 2:
        def getPDBC(x,tag):
            if (x[0] < eps or x[0] > p.L[0] - eps) and (x[1] < eps or x[1] > p.L[1] - eps):
                return np.array([0.0,0.0,0.0])
            elif x[0] < eps or x[0] > p.L[0] - eps:
                return np.array([0.0,round(x[1],5),0.0])
    else:
        def getPDBC(x,flag):
            if (x[0] < eps or x[0] > p.L[0] - eps) and (x[1] < eps or x[1] > p.L[1] - eps) and (x[2] < eps or x[2] > p.L[2] - eps):#x,y,z corner
                return np.array([0.0,0.0,0.0])
            elif (x[0] < eps or x[0] > p.L[0] - eps) and (x[2] < eps or x[2] > p.L[2] - eps):#x-z edge
                return np.array([0.0,round(x[1],5),0.0])
            elif (x[0] < eps or x[0] > p.L[0] - eps) and (x[1] < eps or x[1] > p.L[1] - eps):#x-y edge
                return np.array([0.0,0.0,round(x[2],5)])
            elif (x[1] < eps) or (x[1] > p.L[1]-eps):#on front or back
                return np.array([round(x[0],5),0.0,round(x[2],5)])
            elif (x[0] < eps) or (x[0] > p.L[0]-eps):#on inflow or outflow (left/right)
                return np.array([0.0,round(x[1],5),round(x[2],5)])
    p.periodicDirichletConditions = {0:getPDBC,
                                     1:getPDBC,
                                     2:getPDBC,
                                     3:getPDBC}

pSol = AnalyticalSolutions.PlanePoiseuilleFlow_p(plateSeperation=h,
                                                 mu = mu,
                                                 grad_p = -G)
uSol = AnalyticalSolutions.PlanePoiseuilleFlow_u(plateSeperation=h,
                                                 mu = mu,
                                                 grad_p = -G)
vSol = AnalyticalSolutions.PlanePoiseuilleFlow_v(plateSeperation=h,
                                                 mu = mu,
                                                 grad_p = -G)

p.analyticalSolution = {0:pSol, 1:uSol, 2: vSol}
if p.nd == 3:
    p.analyticalSolution[3] = AnalyticalSolutions.PlanePoiseuilleFlow_v(plateSeperation=h,
                                                                        mu = mu,
                                                                        grad_p = -G)

initialConditions = p.analyticalSolution

nsave=100
dt_init = 1.0e-3
DT = (p.T-dt_init)/float(nsave-1)
p.tnList = [0.0,dt_init]+[dt_init+i*DT for i in range(nsave)]

if  opts.coord:
    if p.nd == 3:
        def onLeft(x):
            return x[0] < eps and x[2] > eps and x[2] < p.L[2] - eps
        def onRight(x):
            return x[0] > p.L[0] - eps and x[2] > eps and x[2] < p.L[2] - eps
        def onFront(x):
            return x[1] < eps and x[2] > eps and x[2] < p.L[2] - eps and x[0] > eps and x[0] < p.L[0] - eps
        def onBack(x):
            return x[1] > p.L[1] - eps and x[2] > eps and x[2] < p.L[2] - eps and x[0] > eps and x[0] < p.L[0] - eps
        def onBottom(x):
            return x[2] < eps
        def onTop(x):
            return x[2] > p.L[2] - eps
    elif p.nd == 2:
        def onLeft(x):
            return x[0] < eps and x[1] > eps and x[1] < p.L[1] - eps
        def onRight(x):
            return x[0] > p.L[0] - eps and x[1] > eps and x[1] < p.L[1] - eps
        def onBottom(x):
            return x[1] < eps
        def onTop(x):
            return x[1] > p.L[1] - eps

if opts.periodic:
    def getDBC_pressure_duct(x,flag):
        pass

    def getDBC_u_duct(x,flag):
        if onTop(x) or onBottom(x):
            return lambda x,t: 0.0

    def getDBC_v_duct(x,flag):
        if onTop(x) or onBottom(x):
            return lambda x,t: 0.0

    def getDBC_w_duct(x,flag):
        if onTop(x) or onBottom(x):
            return lambda x,t: 0.0

    p.dirichletConditions = {0:getDBC_pressure_duct,
                             1:getDBC_u_duct,
                             2:getDBC_v_duct}
    if opts.nd==3:
        p.dirichletConditions[3] = getDBC_w_duct
    

    def getAFBC_p_duct(x,flag):
        if onTop(x) or onBottom(x):
            return lambda x,t: 0.0
        else:
            return lambda x,t: 0.0

    def getAFBC_u_duct(x,flag):
        if onTop(x) or onBottom(x):
            return lambda x,t: 0.0
        else:
            return lambda x,t: 0.0

    def getAFBC_v_duct(x,flag):
        if onTop(x) or onBottom(x):
            return lambda x,t: 0.0
        else:
            return lambda x,t: 0.0

    def getAFBC_w_duct(x,flag):
        if onTop(x) or onBottom(x):
            return lambda x,t: 0.0
        else:
            return lambda x,t: 0.0

    p.advectiveFluxBoundaryConditions =  {0:getAFBC_p_duct,
                                          1:getAFBC_u_duct,
                                          2:getAFBC_v_duct}

    if opts.nd==3:
        p.advectiveFluxBoundaryConditions[3] = getAFBC_w_duct

    def getDFBC_duct(x,flag):
        if onTop(x) or onBottom(x):
            return None
        else:
            return lambda x,t: 0.0

    p.diffusiveFluxBoundaryConditions = {0:{},
                                         1:{1:getDFBC_duct},
                                         2:{2:getDFBC_duct}}

    if opts.nd==3:
        p.diffusiveFluxBoundaryConditions[3] = {3:getDFBC_duct}

else:
    if  opts.coord:
        def getDBC_pressure_duct(x,flag):
            if onRight(x):
                return lambda x,t: 0.0
            
        def getDBC_u_duct(x,flag):
            if onLeft(x):
                return lambda x,t: uSol.uOfX(x)
            if opts.weak and onRight(x):
                return lambda x,t: uSol.uOfX(x)
            if onTop(x) or onBottom(x):
                return lambda x,t: uSol.uOfX(x)
            
        def getDBC_v_duct(x,flag):
            if onLeft(x) or onRight(x) or onTop(x) or onBottom(x):
                return lambda x,t: 0.0

        def getDBC_w_duct(x,flag):
            if onLeft(x) or onRight(x) or onTop(x) or onBottom(x):
                return lambda x,t: 0.0

        p.dirichletConditions = {0:getDBC_pressure_duct,
                                 1:getDBC_u_duct,
                                 2:getDBC_v_duct}
        if p.nd == 3:
            p.dirichletConditions[3] =  getDBC_w_duct


        def getAFBC_p_duct(x,flag):
            if onLeft(x):
                return lambda x,t: -uSol.uOfX(x)
            if onTop(x) or onBottom(x):
                return lambda x,t: 0.0
            if p.nd == 3:
                if onFront(x) or onBack(x):
                    return lambda x,t: 0.0

        def getAFBC_u_duct(x,flag):
            if p.nd == 3:
                if onFront(x) or onBack(x):
                    return lambda x,t: 0.0

        def getAFBC_v_duct(x,flag):
            if p.nd == 3:
                if onFront(x) or onBack(x):
                    return lambda x,t: 0.0

        p.advectiveFluxBoundaryConditions =  {0:getAFBC_p_duct,
                                              1:getAFBC_u_duct,
                                              2:getAFBC_v_duct}
        if p.nd == 3:
            def getAFBC_w_duct(x,flag):
                if onFront(x) or onBack(x):
                    return lambda x,t: 0.0
            p.advectiveFluxBoundaryConditions[3] = getAFBC_w_duct

        def getDFBC_u_duct(x,flag):
            if onRight(x):
                return lambda x,t: 0.0
            elif p.nd == 3:
                if onFront(x) or onBack(x):
                    return lambda x,t: 0.0

        def getDFBC_vw_duct(x,flag):
            if p.nd == 3:
                if onFront(x) or onBack(x):
                    return lambda x,t: 0.0

        p.diffusiveFluxBoundaryConditions = {0:{},
                                             1:{1:getDFBC_u_duct},
                                             2:{2:getDFBC_vw_duct}}
        if p.nd == 3:
            p.diffusiveFluxBoundaryConditions[3] = {3:getDFBC_vw_duct}
    else:
        def getDBC_pressure_duct(x,flag):
            if flag == boundaryTags['right']:
                return lambda x,t: 0.0

        def getDBC_u_duct(x,flag):
            if flag == boundaryTags['left']:
                return lambda x,t: uSol.uOfX(x)
            if opts.weak and flag == boundaryTags['right']:
                return lambda x,t: 0.0
            if flag in [boundaryTags['top'], boundaryTags['bottom']]:
                return lambda x,t: 0.0
        def getDBC_v_duct(x,flag):
            if flag in [boundaryTags['left'],
                        boundaryTags['right'],
                        boundaryTags['top'],
                        boundaryTags['bottom']]:
                return lambda x,t: 0.0

        def getDBC_w_duct(x,flag):
            if flag in [boundaryTags['left'],
                        boundaryTags['right'],
                        boundaryTags['top'],
                        boundaryTags['bottom']]:
                return lambda x,t: 0.0

        p.dirichletConditions = {0:getDBC_pressure_duct,
                                 1:getDBC_u_duct,
                                 2:getDBC_v_duct}
        if p.nd == 3:
            p.dirichletConditions[3] =  getDBC_w_duct


        def getAFBC_p_duct(x,flag):
            if flag == boundaryTags['left']:
                return lambda x,t: -uSol.uOfX(x)
            if flag in [boundaryTags['top'],
                        boundaryTags['bottom']]:
                return lambda x,t: 0.0
            if p.nd == 3 and flag in [boundaryTags['front'],
                                    boundaryTags['back'],
                                    0]:
                return lambda x,t: 0.0
            elif p.nd == 2 and flag == 0:
                return lambda x,t: 0.0

        def getAFBC_u_duct(x,flag):
            if p.nd == 3 and flag in [boundaryTags['front'],
                                    boundaryTags['back'],
                                    0]:
                return lambda x,t: 0.0
            elif p.nd == 2 and flag == 0:
                return lambda x,t: 0.0

        def getAFBC_v_duct(x,flag):
            if p.nd == 3 and flag in [boundaryTags['front'],
                                    boundaryTags['back'],
                                    0]:
                return lambda x,t: 0.0
            elif p.nd == 2 and flag == 0:
                return lambda x,t: 0.0
        def getAFBC_w_duct(x,flag):
            if flag in [boundaryTags['front'],
                        boundaryTags['back'],
                        0]:
                return lambda x,t: 0.0

        p.advectiveFluxBoundaryConditions =  {0:getAFBC_p_duct,
                                            1:getAFBC_u_duct,
                                            2:getAFBC_v_duct}
        if p.nd == 3:
            p.advectiveFluxBoundaryConditions[3] = getAFBC_w_duct

        def getDFBC_duct(x,flag):
            if flag == boundaryTags['right']:#outflow
                return lambda x,t: 0.0
            if p.nd == 3 and flag in [boundaryTags['front'],
                                    boundaryTags['back'],
                                    0]:
                return lambda x,t: 0.0
            elif p.nd == 2 and flag == 0:
                return lambda x,t: 0.0
        p.diffusiveFluxBoundaryConditions = {0:{},
                                             1:{1:getDFBC_duct},
                                             2:{2:getDFBC_duct}}
        if p.nd == 3:
            p.diffusiveFluxBoundaryConditions[3] = {3:getDFBC_duct}

##########
#Numerics#
##########

n = Numerics_base()

#time stepping
n.runCFL = 0.9
if opts.timeOrder == 2:
    n.timeIntegration = TimeIntegration.VBDF
    n.timeOrder = 2
elif opts.timeOrder == 1:
    n.timeIntegration = TimeIntegration.BackwardEuler
    n.timeOrder = 1
elif opts.timeOrder == 0:
    n.timeIntegration = TimeIntegration.NoIntegration

n.stepController  = StepControl.Min_dt_cfl_controller
n.systemStepExact = False

useTaylorHood = False
if opts.spaceOrder == 1:
    if opts.triangles:
        if useTaylorHood:
            n.femSpaces = {0:FemTools.C0_AffineLinearOnSimplexWithNodalBasis,
                           1:FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis,
                           2:FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis}
            if p.nd == 3:
                n.femSpaces[3] = FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis
            n.elementQuadrature = Quadrature.SimplexGaussQuadrature(p.nd,5)
            n.elementBoundaryQuadrature = Quadrature.SimplexGaussQuadrature(p.nd-1,5)
        else:
            n.femSpaces = {0:FemTools.C0_AffineLinearOnSimplexWithNodalBasis,
                           1:FemTools.C0_AffineLinearOnSimplexWithNodalBasis,
                           2:FemTools.C0_AffineLinearOnSimplexWithNodalBasis}
            if p.nd == 3:
                n.femSpaces[3] = FemTools.C0_AffineLinearOnSimplexWithNodalBasis
            n.elementQuadrature = Quadrature.SimplexGaussQuadrature(p.nd,3)
            n.elementBoundaryQuadrature = Quadrature.SimplexGaussQuadrature(p.nd-1,3)
    else:
        if useTaylorHood:
            n.femSpaces = {0:FemTools.C0_AffineLinearOnCubeWithNodalBasis,
                           1:FemTools.C0_AffineQuadraticOnCubeWithNodalBasis,
                           2:FemTools.C0_AffineQuadraticOnCubeWithNodalBasis}
            if p.nd == 3:
                n.hex = True
                n.femSpaces[3] = FemTools.C0_AffineQuadraticOnCubeWithNodalBasis
            else:
                n.quad = True
            n.elementQuadrature = Quadrature.CubeGaussQuadrature(p.nd,2)
            n.elementBoundaryQuadrature = Quadrature.CubeGaussQuadrature(p.nd-1,2)
        else:
            n.femSpaces = {0:FemTools.C0_AffineLinearOnCubeWithNodalBasis,
                           1:FemTools.C0_AffineLinearOnCubeWithNodalBasis,
                           2:FemTools.C0_AffineLinearOnCubeWithNodalBasis}
            if p.nd == 3:
                n.hex = True
                n.femSpaces[3] = FemTools.C0_AffineLinearOnCubeWithNodalBasis
            else:
                n.quad = True
            n.elementQuadrature = Quadrature.CubeGaussQuadrature(p.nd,2)
            n.elementBoundaryQuadrature = Quadrature.CubeGaussQuadrature(p.nd-1,2)

elif opts.spaceOrder == 2:    
    if opts.triangles:
        n.femSpaces = {0:FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis,
                       1:FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis,
                       2:FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis}
        if p.nd == 3:
            n.femSpaces[3] = FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis
        n.elementQuadrature = Quadrature.SimplexGaussQuadrature(p.nd,5)
        n.elementBoundaryQuadrature = Quadrature.SimplexGaussQuadrature(p.nd-1,5)
    else:
        n.femSpaces = {0:FemTools.C0_AffineQuadraticOnCubeWithNodalBasis,
                       1:FemTools.C0_AffineQuadraticOnCubeWithNodalBasis,
                       2:FemTools.C0_AffineQuadraticOnCubeWithNodalBasis}
        if p.nd == 3:
            n.femSpaces[3] = FemTools.C0_AffineQuadraticOnCubeWithNodalBasis
        n.elementQuadrature = Quadrature.SimplexGaussQuadrature(p.nd,5)
        n.elementBoundaryQuadrature = Quadrature.SimplexGaussQuadrature(p.nd-1,5)

n.nnx=opts.nnx
n.nny=(n.nnx-1)//4
if p.nd == 3:
    n.nnz=n.nny

he = p.L[0]/float(n.nnx-1)
if p.nd == 3:
    n.triangleOptions="VApq1.25q12feena%e" % ((he**3)/6.0,)
else:
    n.triangleFlag = 1#if regular triangulatio then alternate diagonals
    n.triangleOptions="pAq30.0Dena%f" % ((he**2)/4.0,)

n.numericalFluxType = RANS2P.NumericalFlux

if opts.periodic:
    n.periodicDirichletConditions = p.periodicDirichletConditions
    n.parallelPeriodic=True

n.subgridError = RANS2P.SubgridError(coefficients=p.coefficients,
                                     nd=p.nd,
                                     lag=False,
                                     hFactor=1.0)

n.shockCapturing = RANS2P.ShockCapturing(coefficients=p.coefficients,
                                         nd=p.nd,
                                         shockCapturingFactor=0.0,
                                         lag=True)

n.multilevelNonlinearSolver  = NonlinearSolvers.Newton

n.levelNonlinearSolver = NonlinearSolvers.Newton

n.fullNewtonFlag = True
n.maxNonlinearIts = 50
n.maxLineSearches = 0

n.tolFac = 0.0

n.nl_atol_res = 1.0e-8

n.matrix = LinearAlgebraTools.SparseMatrix

n.multilevelLinearSolver = LinearSolvers.KSP_petsc4py
n.levelLinearSolver = LinearSolvers.KSP_petsc4py

if opts.pc_type == 'LU':
    n.multilevelLinearSolver = LinearSolvers.LU
    n.levelLinearSolver = LinearSolvers.LU
    n.linearSmoother = None
elif opts.pc_type == 'selfp_petsc':
    if p.nd==3:
        n.linearSmoother = LinearSolvers.SimpleNavierStokes3D
        n.linearSmootherOptions = (opts.A_block_AMG,)
    elif p.nd==2:
        n.linearSmoother = LinearSolvers.SimpleNavierStokes2D
        n.linearSmootherOptions = (opts.A_block_AMG,)
elif opts.pc_type == 'two_phase_PCD':
    n.linearSmoother = LinearSolvers.NavierStokes_TwoPhasePCD
    n.linearSmootherOptions = (False,
                               True,
                               1,
                               0,
                               1,
                               False)
    #(density_scaling, numerical_viscosity, pcd_lumped, chebyshev_its, laplace_null_space)

n.linear_solver_options_prefix = 'rans2p_'

n.linTolFac = 0.1
n.l_atol_res = 0.1*n.nl_atol_res

n.conservativeFlux = None

if opts.periodic:
    n.parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.element
else:
    n.parallelPartitioningType = MeshTools.MeshParallelPartitioningTypes.node
n.nLayersOfOverlapForParallel = 0

##############
#System input#
##############
from proteus.default_so import *
systemStepExact=n.systemStepExact
tnList = p.tnList
pnList=[(p,n)]
if opts.periodic:
    mesh_name = "pg"
else:
    if opts.grid:
        mesh_name = "rg"
    else:
        mesh_name = "ug"
if opts.triangles:
    space_name="p{0}".format(opts.spaceOrder)
else:
    space_name="1{0}".format(opts.spaceOrder)

name = "duct{0}t{1}{2}d{3}he{4}".format(space_name,
                                               opts.timeOrder,
                                               opts.nd,
                                               mesh_name,
                                               he)
