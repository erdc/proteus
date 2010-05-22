#!/usr/bin/env python
from EGeometry import *
from FemTools import *
from MeshTools import *
from LinearAlgebraTools import *
from LinearSolvers import *
from NonlinearSolvers import *
from QuadTools import *
from LevelSetTests import *
from TimeIntegrationTools import *
import ScalarTransport
import AnalyticalSolutions
## \defgroup RKDG RKDG
#
# @{

"""
define some of the finite element tools for RKDG discretizations

TODO:

"""
class ShockIC:
    def uOfX(self,X):
        if X <= 0.5:
            return 1.0
        else:
            return 0.0
    def uOfXT(self,X,T):
        if X <= 0.5:
            return 1.0
        else:
            return 0.0
#end shockIC
def getDBC_shock(x):
    if x == 0.0:
        return lambda x,t: 1.0
    if x == 1.0:
        return lambda x,t: 0.0
#end def
class SlugIC:
    def uOfX(self,X):
        if 0.1 <= X and X <= 0.4:
            return 1.0
        else:
            return 0.0
    def uOfXT(self,X,T):
        if 0.1 <= X and X <= 0.4:
            return 1.0
        else:
            return 0.0
#end shockIC


def testAdv1d(mesh,
              runCFL=0.1,
              ntMax=1000,
              femOrder=1,
              tIntOrder=2,
              initCond=None,
              bndyCond=None,
              quadratureOrder=3,
              bquadratureOrder=3,
              errorQuadratureOrder=3,
              verbose=0):
    """
    run 1d dirac delta problem with RKDG 

    """
    #mwf debug
    print 'entering testAdv1d args are'
    print 'runCFL = ',runCFL
    print 'ntMax= ',ntMax
    print 'femOrder= ',femOrder
    print 'tIntOrder= ',tIntOrder
    print 'initCond= ',initCond
    print 'bndyCond= ',bndyCond
    print 'quadratureOrder= ',quadratureOrder
    print 'bquadratureOrder= ',quadratureOrder
    print 'errorQuadratureOrder= ',quadratureOrder
    print 'verbose= ',verbose
    
    nd = 1
    T = 0.3
    t0= 0.0
    ntFix = 10
    if verbose > 0:
        print 'definining problem coefficients'
    #end verbose
    A0=Numeric.array([[0.0]])
    B0=Numeric.array([1.0])
    C0=0.0

    exactSoln    = AnalyticalSolutions.LinearADR_Decay_DiracIC(b=B0,a=1.0e-3,c=C0,tStart=0.25)
    probCoefficients = ScalarTransport.LinearADR_ConstantCoefficients(M=1.0,A=A0,B=B0,C=C0)
    probCoefficients.mass = 'linear'
    probCoefficients.advection = 'linear'
    probCoefficients.diffusion = None
    probCoefficients.potential = None
    probCoefficients.reaction  = None
    #probCoefficients.diffusion = 'constant'
    #probCoefficients.potential = 'linear'
    #probCoefficients.reaction  = 'linear'
    
    if bndyCond == None:
        bndyCond  = getDBC_hom
    if initCond == None:
        initCond  = exactSoln
    #end if
    fluxBndyCond = 'outFlow'
    
    fullNewtonFlag= False
    maxNLits  = 2
    
    if verbose > 0:
        print 'definining quadrature'
    #end verbose
    quadrature = {}
    elemQuadrature = SimplexGaussQuadrature(nd)
    elemQuadrature.setOrder(quadratureOrder)
    for integral in ScalarTransport.OneLevelScalarTransport.integralKeys:
        quadrature[integral] = elemQuadrature
    if useCG:
        quadrature['stab']=elemQuadrature
        quadrature['numDiff']=elemQuadrature
    #end for
    elementBoundaryQuadrature = {}
    boundaryQuadrature = SimplexGaussQuadrature(nd-1)
    boundaryQuadrature.setOrder(bquadratureOrder)
    for bintegral in ScalarTransport.OneLevelScalarTransport.elementBoundaryIntegralKeys:
        elementBoundaryQuadrature[bintegral] = boundaryQuadrature
    # # setup finite element spaces and stuff
    if verbose > 0:
        print 'definining finite element spaces'
    #end verbose
    #try P^2 too
    print 'rotating cone using dg order= ',femOrder
    print 'rotating cone using rk order= ',tIntOrder
    if femOrder == 1:
        FemSpace = DG_AffineLinearOnSimplexWithNodalBasis(mesh,nd)
    else:
        FemSpace = DG_AffineQuadraticOnSimplexWithNodalBasis(mesh,nd)
    #end if
    u   = FiniteElementFunction(FemSpace)
    phi = FiniteElementFunction(FemSpace)

    if verbose > 0:
        print 'definining boundary conditions'
    #end verbose
    dirichletBCs=DOFBoundaryConditions(FemSpace,bndyCond)

    # # setup numerical approximation configuration
    if verbose > 0:
        print 'definining linear algebra'
    #end verbose
    MatType = SparseMat
    matType = 'csr'
    linearSolverType=  'SparseLU'
    # # time integration
    #steady state
    if verbose > 0:
        print 'definining time integration'
    #end verbose
    TimeIntegrationClass = SSPRKintegration
    #mwf debug
    #TimeIntegrationClass = BackwardEuler
    # # flux approximations
    if verbose > 0:
        print 'definining flux approximations'
    #end verbose
    
    #pick fluxes for DG
    numericalFlux = True
    conservativeFlux = None
    stabilization = None
    shockCapturing=None
    shockCapturingDiffusion=None
    
    # #create a single level solver
    if verbose > 0:
        print 'creating single level system'
    #end verbose
    system = ScalarTransport.OneLevelScalarTransport(u,
                                                     phi,
                                                     FemSpace,
                                                     dirichletBCs,
                                                     probCoefficients,
                                                     quadrature,
                                                     elementBoundaryQuadrature,
                                                     fluxBndyCond,
                                                     stabilization,
                                                     shockCapturing,
                                                     shockCapturingDiffusion,                
                                                     conservativeFlux,
                                                     numericalFlux,
                                                     TimeIntegrationClass,
                                                     tIntOrder)

    #need this?
    if verbose > 0:
        print 'femSpace dim= ',FemSpace.dim
        print 'system dim = ',system.dim
    #end if
    y  = Numeric.zeros((system.dim,),Numeric.Float)
    dy = Numeric.zeros((system.dim,),Numeric.Float)
    r  = Numeric.zeros((system.dim,),Numeric.Float)
    system.timeIntegration.chooseDT()
    system.updateQuadrature()
    system.setFreeDOF(y)
    system.updateCoefficients()
    system.updateQuadrature()
    system.setInitialConditions(initCond,t0)
    system.setFreeDOF(y)
    if verbose > 0:
        print 'after setInitialConditions'
        print 'u.dof= ',u.dof
    #end if
    #mwf debug
    #print """before system.timeIntegration.chooseDT runCFL= %g DT = %g
    #timeIntegration.cfl = %s
    #""" % (system.timeIntegration.runCFL,
    #       system.timeIntegration.DT,
    #       system.timeIntegration.cfl)
    
    system.timeIntegration.chooseDT()
    
    #mwf debug
    #print """before getResidual system.timeIntegration.DT = %g
    #""" % system.timeIntegration.DT
    system.getResidual(y,r)
    
    if verbose > 0:
        print 'y0= ',y
        print 'initial residual = ',r
        print 'initial dt= ',system.timeIntegration.DT
        print 'u.dof= ',u.dof
    #end verbose
    if nd == 1:
        FemSpace.writeFunctionGnuplot(u,'init')
        raw_input('Please press return to continue... \n')
        
    else:
        FemSpace.writeFunctionMatlab(u,'init','mesh2drkdg')
    
    # # create nonlinear system for solving problem
    if verbose > 0:
        print 'creating nonlinear system'
    #end verbose
    
    #jacobian = MatType(system.dim,system.dim,min(7,system.dim))
    nblock = u.femSpace.referenceFiniteElement.localFunctionSpace.dim
    jacobian = MatType(system.dim,system.dim,min(nblock,system.dim))
    #cek set up sparse  matrix  structures
    system.matType = matType
    jacobian = system.initializeJacobian(jacobian)
    # # create linear system for solving problem
    if verbose > 5:
        print 'creating linear system'
        print 'initial jacobian mat = ',jacobian
    #end verbose
    
    linearSolver = SparseLU(jacobian)
    
    nlSolver = Newton(linearSolver,system,jacobian,maxIts=maxNLits,
                      fullNewton=fullNewtonFlag,printInfo=verbose > 2)
    #need intial residual calculated
    
    system.timeIntegration.runCFL = runCFL
    t = t0
    nsteps= 0
    nstages = tIntOrder
    dtFix = float(T/ntFix)
    adaptDT = False
    while t < T and nsteps < ntMax:
        
        system.timeIntegration.chooseDT()
        dtMin = min(T-t,system.timeIntegration.DT)
        if not adaptDT:
            dtMin = min(dtMin,dtFix)
        #end not adapting
        if T-(t+dtMin) <= 1.0e-10:
            dtMin = T-t
        #end dt set
        if dtMin <= 1.0e-14:
            print 'dtMin = ',dtMin, ' too small quitting '
            sys.exit(1)
        #end if
        system.timeIntegration.DT=dtMin
        if nsteps == 0:
            system.initializeTimeIntegration()
        else:
            system.timeIntegration.updateTimeHistory()
        #end if
        t += system.timeIntegration.DT #do this before or after solve?
        
        for i in range(nstages):
            system.getResidual(y,r)
            #solve system for new time level
            nlSolver.solve(y,r)
            #have to have this here for multistage values?
            #system.updateCoefficients()
            system.timeIntegration.updateStage()
            
        #end stage loop
        nsteps += 1

        #mwf prevent too small steps?
        #if t >= T-1.0e-12:
        #    t = T
        
        if verbose > -1:
            print 'took step to t= ',t,' DT= ',system.timeIntegration.DT
        if nd == 1 and verbose > -1:
            FemSpace.writeFunctionGnuplot(u,'solution')
        if verbose > 5:
            nx = mesh.nElements_global
            ndof = u.femSpace.referenceFiniteElement.localFunctionSpace.dim
            print 'u.dof= ',Numeric.reshape(u.dof,(nx,ndof))
        #end if
        #end if verbose
        #make sure solution is kept by ode
        system.setUnknowns(y)
        
        if verbose > 5:
            print 'u.dof= ',u.dof
            print 'y= ',y
    #end time loop
    #go ahead and get coefficients at final time solution
    system.timeIntegration.updateTimeHistory()
    system.updateCoefficients()
    print 'reached t= ',t,' nsteps= ',nsteps
    if nd == 1:
        FemSpace.writeFunctionGnuplot(u,'solution')
    else:
        FemSpace.writeFunctionMatlab(u,'solution','mesh2drkdg')
    #end if
    raw_input('Please press return to continue... \n')

def testRotatingCone2d(mesh,runCFL=0.1,
                       femOrder=1,tIntOrder=2,
                       quadratureOrder=4,
                       bquadratureOrder=3,
                       errorQuadratureOrder=3,
                       ntMax=1000,
                       verbose=0):
    """
    run 2d rotating cone problem with RKDG 

    """
    #mwf debug
    print 'entering testRotatingCone2d args are'
    print 'runCFL = ',runCFL
    print 'ntMax= ',ntMax
    print 'femOrder= ',femOrder
    print 'tIntOrder= ',tIntOrder
    print 'quadratureOrder= ',quadratureOrder
    print 'bquadratureOrder= ',bquadratureOrder
    print 'errorQuadratureOrder= ',errorQuadratureOrder
    print 'verbose= ',verbose

    nd = 2
    T = 0.5
    t0= 0.0
    ntFix = 10
    if verbose > 0:
        print 'definining problem coefficients'
    #end verbose
    radius = 0.125
    exactSoln    = RotatingCone2D(radius)
    probCoefficients = UnitSquareRotation()
    probCoefficients.mass = 'linear'
    probCoefficients.advection = 'linear'
    probCoefficients.diffusion = None
    probCoefficients.potential = None
    probCoefficients.reaction  = None
    
#    bndyCond  = getHomogeneousDBC2D
    bndyCond  = getNoDBC
    fluxBndyCond = 'outFlow'
    initCond  = exactSoln
    fullNewtonFlag= False
    maxNLits  = 1
    
    if verbose > 0:
        print 'definining quadrature'
    #end verbose
    quadrature = {}
    elemQuadrature = SimplexGaussQuadrature(nd)
    elemQuadrature.setOrder(quadratureOrder)
    for integral in ScalarTransport.OneLevelScalarTransport.integralKeys:
        quadrature[integral] = elemQuadrature
    if useCG:
        quadrature['stab']=elemQuadrature
        quadrature['numDiff']=elemQuadrature
    #end for
    elementBoundaryQuadrature = {}
    boundaryQuadrature = SimplexGaussQuadrature(nd-1)
    boundaryQuadrature.setOrder(bquadratureOrder)
    for bintegral in ScalarTransport.OneLevelScalarTransport.elementBoundaryIntegralKeys:
        elementBoundaryQuadrature[bintegral] = boundaryQuadrature
    # # setup finite element spaces and stuff
    if verbose > 0:
        print 'definining finite element spaces'
    #end verbose
    #try P^2 too
    print 'rotating cone using dg order= ',femOrder
    print 'rotating cone using rk order= ',tIntOrder
    if femOrder == 1:
        FemSpace = DG_AffineLinearOnSimplexWithNodalBasis(mesh,nd)
    else:
        FemSpace = DG_AffineQuadraticOnSimplexWithNodalBasis(mesh,nd)
    #end if
    u   = FiniteElementFunction(FemSpace)
    phi = FiniteElementFunction(FemSpace)
    ndofLoc = FemSpace.referenceFiniteElement.localFunctionSpace.dim
    if verbose > -1:
        print 'local function space dim= ',ndofLoc
    #end verbose
    
    if verbose > 0:
        print 'definining boundary conditions'
    #end verbose
    dirichletBCs=DOFBoundaryConditions(FemSpace,bndyCond)
    
    # # setup numerical approximation configuration
    if verbose > 0:
        print 'definining linear algebra'
    #end verbose
    MatType = SparseMat
    matType = 'csr'
    linearSolverType=  'SparseLU'
    # # time integration
    #steady state
    if verbose > 0:
        print 'definining time integration'
    #end verbose
    TimeIntegrationClass = SSPRKintegration
    #mwf debug
    #TimeIntegrationClass = BackwardEuler
    # # flux approximations
    if verbose > 0:
        print 'definining flux approximations'
    #end verbose
    
    #pick fluxes for DG
    conservativeFlux = None
    numericalFlux = True
    if useCG:
        stabilization = '2'
        shockCapturing=None
        shockCapturingDiffusion=0.1
    else:
        stabilization= None
        shockCapturing=None
        shockCapturingDiffusion=None

    # #create a single level solver
    if verbose > 0:
        print 'creating single level system'
    #end verbose
    system = ScalarTransport.OneLevelScalarTransport(u,
                                                     phi,
                                                     FemSpace,
                                                     dirichletBCs,
                                                     probCoefficients,
                                                     quadrature,
                                                     elementBoundaryQuadrature,
                                                     fluxBndyCond,
                                                     stabilization,
                                                     shockCapturing,
                                                     shockCapturingDiffusion,                
                                                     conservativeFlux,
                                                     numericalFlux,
                                                     TimeIntegrationClass,
                                                     tIntOrder)

    #need this?
    if verbose > 0:
        print 'femSpace dim= ',FemSpace.dim
        print 'system dim = ',system.dim
    #end if
    y  = Numeric.zeros((system.dim,),Numeric.Float)
    dy = Numeric.zeros((system.dim,),Numeric.Float)
    r  = Numeric.zeros((system.dim,),Numeric.Float)
    system.timeIntegration.chooseDT()
    system.updateQuadrature()
    system.setFreeDOF(y)
    system.updateCoefficients()
    system.updateQuadrature()
    system.setInitialConditions(initCond,t0)
    system.setFreeDOF(y)
    if verbose > 0:
        print 'after setInitialConditions'
        print 'u.dof= ',u.dof
    #end if
    #mwf debug
    #print """before system.timeIntegration.chooseDT runCFL= %g DT = %g
    #timeIntegration.cfl = %s
    #""" % (system.timeIntegration.runCFL,
    #       system.timeIntegration.DT,
    #       system.timeIntegration.cfl)
    system.timeIntegration.chooseDT()

    system.getResidual(y,r)
    
    if verbose > 0:
        print 'y0= ',y
        print 'initial residual = ',r
        print 'initial dt= ',system.timeIntegration.DT
        print 'u.dof= ',u.dof
    #end verbose
    FemSpace.writeFunctionMatlab(u,'init','mesh2drkdg')

    # # create nonlinear system for solving problem
    if verbose > 0:
        print 'creating nonlinear system'
    #end verbose
    jacobian = MatType(system.dim,system.dim,min(ndofLoc,system.dim))
    #cek set up sparse  matrix  structures
    system.matType = matType
    jacobian = system.initializeJacobian(jacobian)
    # # create linear system for solving problem
    if verbose > 5:
        print 'creating linear system'
        print 'initial jacobian mat = ',jacobian
        print 'maxNLits=',maxNLits
    #end verbose

    linearSolver = SparseLU(jacobian)

    nlSolver = Newton(linearSolver,system,jacobian,maxIts=maxNLits,
                      fullNewton=fullNewtonFlag,printInfo=verbose > 1)


    system.timeIntegration.runCFL = runCFL
    t = t0
    nsteps= 0
    nstages = tIntOrder
    dtFix = float(T/ntFix)
    adaptDT = False
    while t < T and nsteps < ntMax:
        
        system.timeIntegration.chooseDT()
        dtMin = min(T-t,system.timeIntegration.DT)
        if not adaptDT:
            dtMin = min(dtMin,dtFix)
        #end not adapting
        if T-(t+dtMin) <= 1.0e-10:
            dtMin = T-t
        #end dt set
        if dtMin <= 1.0e-14:
            print 'dtMin = ',dtMin, ' too small quitting '
            sys.exit(1)
        #end if

        system.timeIntegration.DT=dtMin
        if nsteps == 0:
            if verbose > -1:
                print 'initailizing on first step'
            #end if
            system.initializeTimeIntegration()
        else:
            system.timeIntegration.updateTimeHistory()
        #end if
        t += system.timeIntegration.DT #do this before or after solve?
        
        for i in range(nstages):
            if verbose > -1:
                print 'step to ',t,' stage solve number ',i
            #end if
            if verbose > 5:
                print 'jacobian= ',nlSolver.J
            #end if
            system.getResidual(y,r)
            #solve system for new time level
            nlSolver.solve(y,r)
            #have to have this here for multistage values?
            #system.updateCoefficients()
            system.timeIntegration.updateStage()
            if verbose > 1:
                print 'after stage ',i
                print 'y= ',y
            #end if
        #end stage loop
        nsteps += 1

        if verbose > -1:
            print 'took step to t= ',t,' DT= ',system.timeIntegration.DT
        if verbose > 0:
            FemSpace.writeFunctionGnuplot(u,'solution')

        #end if verbose
        #make sure solution is kept by ode
        system.setUnknowns(y)

        if verbose > 5:
            print 'u.dof= ',u.dof
            print 'y= ',y
    #end time loop
    #go ahead and get coefficients at final time solution
    system.timeIntegration.updateTimeHistory()
    system.updateCoefficients()
    print 'reached t= ',t,' nsteps= ',nsteps
    if verbose > 2:
        print 'u.dof= ',u.dof
    #end if
    FemSpace.writeFunctionMatlab(u,'solution','mesh2drkdg')

    raw_input('Please press return to continue... \n')
    
## @}

if __name__ == '__main__':

    import sys
    from MeshTools import TriangularMesh
    from PoissonTestProblems import *
    import DiagUtils
    
    from optparse import OptionParser
    parser = OptionParser()
    #options controlling simulation behavior
    parser.add_option('-d','--spacedim',
                      default=1,
                      help="""number of spatial dimensions [1]""")
    parser.add_option('--useDG',
                      default=False,action='store_true',
                      help="""use DG approximation [False]""")
    parser.add_option('--ntMax',
                      default=1000,
                      help="""max number of time steps allowed [1000]""")
    parser.add_option('--nx',
                      default=5,
                      help="""number of nodes along x axis [5]""")
    parser.add_option('--ny',
                      default=5,
                      help="""number of nodes along y axis [5]""")

    parser.add_option('--nz',
                      default=5,
                      help="""number of nodes along z axis [5]""")

    parser.add_option('-O','--order',
                      default=2,
                      help="""order of spatial approximation to use [2]""")

    parser.add_option('--torder',
                      default=-1,
                      help="""order of temporal approximation to use [order+1]""")

    parser.add_option('--runCFL',
                      default=0.1,
                      help="""target CFL for transient simulation [0.1]""")


    parser.add_option('-t','--testNumber',
                      default=0,
                      help="""which test to use [0]
                      0 --- testBarycentricCoords
                      1 --- testQuadNodalBasis
                      2 --- testQuadDOFMap
                      3 --- testLaplacian
                      4 --- rotating cone2d
                      5 --- 1d adv: dirac delta
                      6 --- 1d adv: Riemann ic
                      7 --- 1d adv: slug ic
                      """)
    parser.add_option('--viewMesh',
                      default=False,action='store_true',
                      help="""look at mesh [False]""")

    parser.add_option('-v','--verbose',
                      default=0,
                      help="""level of verbosity in simulation [0]""")

    #get options
    options, args = parser.parse_args(sys.argv[1:]) #get command line args
    #
    verbose = int(options.verbose)
    useDG   = options.useDG
    viewMesh= options.viewMesh
    order   = int(options.order)
    tIntOrder = int(options.torder)
    #mwf debug
    useCG = not useDG
    if verbose > 0:
        print 'main useCG= ',useCG
    #end verbose
    #number of space dimensions
    nd = int(options.spacedim)
    nx = int(options.nx)
    ny = int(options.ny)
    nz = int(options.nz)
    ntMax  = int(options.ntMax)
    runCFL = float(options.runCFL)
    testParams = {}
    testParams['number'] = int(options.testNumber)
    testParams['nd']     = nd
    #need mesh now
    #create mesh
    Lx = 1.0
    Ly = 1.0
    Lz = 1.0
    if nd == 1:
        grid1d = RectangularGrid(nx,1,1,Lx,Ly,Lz)
        mesh1d = EdgeMesh()
        mesh1d.rectangularToEdge(grid1d)
        mesh1d.computeGeometricInfo()
        testParams['mesh'] = mesh1d
        if viewMesh == True:
            testParams['mesh'].writeEdgesGnuplot('mesh1d')
            testParams['mesh'].viewMeshGnuplotPipe('mesh1d')
        #end 
    elif nd == 2:
        mesh2d = TriangularMesh()
        mesh2d.constructTriangularMeshOnRectangle(Lx,Ly,nx,ny)
        mesh2d.computeGeometricInfo()
        testParams['mesh'] = mesh2d
        if viewMesh == True:
            testParams['mesh'].writeEdgesGnuplot('mesh2d')
            testParams['mesh'].viewMeshGnuplotPipe('mesh2d')
        #end 
    elif nd == 3:
        grid3d = RectangularGrid(nx,ny,nz,Lx,Ly,Lz)
        mesh3d = TetrahedralMesh()
        mesh3d.rectangularToTetrahedral(grid3d)
        mesh3d.computeGeometricInfo()
        testParams['mesh'] = mesh3d
        if viewMesh == True:
            testParams['mesh'].writeEdgesGnuplot('mesh3d')
            testParams['mesh'].viewMeshGnuplotPipe('mesh3d')
        #end if
    #end if on nd
    #pick finite element space
    if useCG:
        if order == 1:
            testParams['FemSpace'] = C0_AffineLinearOnSimplexWithNodalBasis(testParams['mesh'],
                                                                            testParams['nd'])
        else:
            testParams['FemSpace'] = C0_AffineQuadraticOnSimplexWithNodalBasis(testParams['mesh'],
                                                                               testParams['nd'])
    else:
        if order == 1:
            testParams['FemSpace'] = DG_AffineLinearOnSimplexWithNodalBasis(testParams['mesh'],
                                                                            testParams['nd'])
        else:
            testParams['FemSpace'] = DG_AffineQuadraticOnSimplexWithNodalBasis(testParams['mesh'],
                                                                               testParams['nd'])
        #end if on order
    #end if on CG for FemSpace

    if testParams['number'] == 0:
        RefUtils.testBarycentricCoords(verbose)
    elif testParams['number'] == 1:
        DiagUtils.testQuadNodalBasis(nd,verbose)
    elif testParams['number'] == 2:
        DiagUtils.testQuadDOFMap(testParams['mesh'],nd,verbose)
    elif testParams['number'] == 3:
        #pick test problem
        if nd == 1:
            poisProb = collectPoissonDescription(uexP5,duexP5,qexP5,AP5,
                                                 DiagUtils.Ident1,fP5,dirBCP5)
        elif nd == 2:
            poisProb = collectPoissonDescription(uexP1,duexP1,qexP1,AP1,
                                                 DiagUtils.Ident2,fP1,dirBCP1)
        else:
            poisProb = collectPoissonDescription(uexP4,duexP4,qexP4,AP4,
                                                 DiagUtils.Ident3,fP4,dirBCP4)
        #end if on dim for poisProb
        quadOrder = 3
        bquadOrder= 3
        equadOrder= 3
        testLaplacian(testParams['mesh'],testParams['FemSpace'],poisProb,
                      Lx,Ly,Lz,nx,ny,nz,nd,useCG,order,
                      quadOrder,bquadOrder,equadOrder,
                      verbose)
    elif testParams['number'] == 4:
        if not nd == 2:
            print 'test number 4, rotating cone in 2d requires nd=2'
            sys.exit(1)
        #end if
        if tIntOrder <= 0:
            tIntOrder = order+1
        #end if
        testRotatingCone2d(mesh2d,runCFL=runCFL,
                           femOrder=order,tIntOrder=tIntOrder,
                           ntMax=ntMax,verbose=verbose)
    elif testParams['number'] == 5:
        if not nd == 1:
            print 'test number 5, linear AD with dirac is requires nd=1'
            sys.exit(1)
        #end if
        if tIntOrder <= 0:
            tIntOrder = order+1
        #end if
        print 'verbose = ',verbose
        testAdv1d(mesh1d,runCFL=runCFL,
                  femOrder=order,tIntOrder=tIntOrder,
                  ntMax=ntMax,verbose=verbose)
    elif testParams['number'] == 6:
        if not nd == 1:
            print 'test number 6, linear AD with heaviside ic requires nd=1'
            sys.exit(1)
        #end if
        if tIntOrder <= 0:
            tIntOrder = order+1
        #end if
        initCond = ShockIC()
        bndyCond = getDBC_shock
        testAdv1d(mesh1d,initCond=initCond,bndyCond=bndyCond,
                  runCFL=runCFL,
                  femOrder=order,tIntOrder=tIntOrder,
                  ntMax=ntMax,verbose=verbose)
    elif testParams['number'] == 7:
        if not nd == 1:
            print 'test number 6, linear AD with slug ic requires nd=1'
            sys.exit(1)
        #end if
        if tIntOrder <= 0:
            tIntOrder = order+1
        #end if
        initCond = SlugIC()
        bndyCond = getDBC_hom
        testAdv1d(mesh1d,initCond=initCond,bndyCond=bndyCond,
                  runCFL=runCFL,
                  femOrder=order,tIntOrder=tIntOrder,
                  ntMax=ntMax,verbose=verbose)
    else:
        print 'test= ',testParams['number'],' not implemented!'
    #end else
