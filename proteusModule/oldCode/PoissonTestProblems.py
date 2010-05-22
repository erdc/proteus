#!/usr/bin/env/python

from EGeometry import *
from ScalarTransport import *
from LinearAlgebraTools import *
from LinearSolvers import *
from NonlinearSolvers import *
from QuadTools import *
import DiagUtils 
import RefUtils 

## \defgroup PoissonTestProblems PoissonTestProblems
#
# @{

class PoissonCoefficients(ScalarTransportCoefficients):
    """
    This class implements:
    -\deld(A\grad u) -S(x) = 0

    """
    def __init__(self,A,S,nd):
        ScalarTransportCoefficients.__init__(self,
                                             mass=None,
                                             advection=None,
                                             diffusion='constant',
                                             potential='linear',
                                             reaction='linear')
        self.A = A
        self.S = S
        self.nSpaceDim = nd
    def evaluate(self,
                 t,
                 x,
                 u,
                 m,dm,
                 f,df,
                 a,da,
                 phi,dphi,
                 r,dr):
        m[:] = 0.0
        dm[:] = 0.0
        
        f[:]=0.0
        df[:] = 0.0
        
        a[:,:] = self.A
        da[:,:] = 0.0

        phi[:] = u[:]
        dphi[:]= 1.0
        for i in range(Numeric.size(x,0)):
            r[i] = -self.S(x[i])
        dr[:]= 0.0

    #end evaluate
#end poisson

class poissonExactSoln:
    """
    just hold solution u(x,t) in format code is expecting
    uses default values if 
    """
    def __init__(self,fxt):
        self.fxt = fxt
    #end init
    def uOfX(self,X):
        return fxt(X,0.0)
    def uOfXT(self,X,T):
        return self.fxt(X,T)
#end poissonExactSoln
def testLaplacian(mesh,FemSpace,poissonProb,
                  Lx,Ly,Lz,nx,ny,nz,
                  nd=2,useCG=False,order=2,
                  quadratureOrder=3,
                  bquadratureOrder=3,
                  errorQuadratureOrder=3,
                  verbose=0):
    """
    try to write some code to define a simple linear steady-state problem and solve it
    using the ScalarTransport interface and CrR elements

    solve
       -\deld u = f   in \Omega = [0,1]^d
              u = u^b   on \partial \Omega
              u = \sum^d_i=1 x_i^2

    """
    import P1ncMixed
    import ScalarTransport
    import TimeIntegrationTools
    # #  setup problem coefficients 
    if verbose > 0:
        print 'definining problem coefficients'
    #end verbose
    
    fProb   = poissonProb['f']
    AProb   = poissonProb['A']
    A0Prob  = poissonProb['A0']
    uexProb = poissonProb['uex']
    qexProb = poissonProb['qex']
    dirProb = poissonProb['DirichletBCs']

    #end else on space dim
    probCoefficients = PoissonCoefficients(A=A0Prob,S=fProb,nd=nd)
    exactSoln = poissonExactSoln(uexProb)
    qexactSoln= poissonExactSoln(qexProb)
    # # setup boundary conditions
    # # setup numerical quadrature
    if verbose > 0:
        print 'definining quadrature'
    #end verbose
    quadrature = {}
    #end quadrature order
    elemQuadrature = SimplexGaussQuadrature(nd)
    elemQuadrature.setOrder(quadratureOrder)
    for integral in ScalarTransport.OneLevelScalarTransport.integralKeys:
        quadrature[integral] = elemQuadrature
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
    #try P^2
    if useCG:
        if order == 1:
            FemSpace = C0_AffineLinearOnSimplexWithNodalBasis(mesh,nd)
        else:
            FemSpace = C0_AffineQuadraticOnSimplexWithNodalBasis(mesh,nd)
    else:
        if order == 1:
            FemSpace = DG_AffineLinearOnSimplexWithNodalBasis(mesh,nd)
        else:
            FemSpace = DG_AffineQuadraticOnSimplexWithNodalBasis(mesh,nd)

    #end if
    u   = FiniteElementFunction(FemSpace)
    phi = FiniteElementFunction(FemSpace)
    #phi = u
    if verbose > 0:
        print 'definining boundary conditions'
    #end verbose
    dirichletBCs=DOFBoundaryConditions(
        FemSpace,dirProb)
    fluxBndyCond='noFlow'
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
    TimeIntegrationClass = TimeIntegrationTools.TimeIntegration
    # # flux approximations
    if verbose > 0:
        print 'definining flux approximations'
    #end verbose
    #diffusion only
    if useCG:
        conservativeFlux = None
        numericalFlux    = None
        stabilization = None
        shockCapturing= None
        shockCapturingDiffusion= None
    else:
        conservativeFlux = None
        numericalFlux    = True
        stabilization = None
        shockCapturing= None
        shockCapturingDiffusion= None
    #end
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
                                                     TimeIntegrationClass)
    #need this?
    if verbose > 0:
        print 'femSpace dim= ',FemSpace.dim
        print 'system dim = ',system.dim
    #end if
    y  = Numeric.zeros((system.dim,),Numeric.Float)
    dy = Numeric.zeros((system.dim,),Numeric.Float)
    r  = Numeric.zeros((system.dim,),Numeric.Float)
    system.updateQuadrature()
    system.setFreeDOF(y)
    system.updateCoefficients()
    system.updateQuadrature()
    system.matType = matType
    
    # # create nonlinear system for solving problem
    if verbose > 0:
        print 'creating nonlinear system'
    #end verbose
    
    jacobian = MatType(system.dim,system.dim,min(7,system.dim))
    #mwf add following Chris example
    jacobian = system.initializeJacobian(jacobian)
    # # create linear system for solving problem
    if verbose > 5:
        print 'creating linear system'
        print 'initial jacobian mat = ',jacobian
    #end verbose

    linearSolver = SparseLU(jacobian)

    nlSolver = Newton(linearSolver,system,jacobian,maxIts=10)
    #need intial residual calculated
    system.getResidual(y,r)
    if verbose > 5:
        print 'system y0= ',y
        system.getJacobian(jacobian)
        print 'system residual= ',r
        print 'system jacobian = ',jacobian
    # # solving problem
    if verbose > 0:
        print 'trying to solve nonlinear system'
    #end verbose
    nlSolver.solve(y,r)

        
    if verbose > 0:
        print 'nonlinear solver done'
    #end if

    #create some numerical quadrature points to evaluate errors. Allow for
    #these to be different order than the ones used in solution method

    errorQuadrature = SimplexGaussQuadrature(nd)
    errorQuadrature.setOrder(errorQuadratureOrder)

    errQuadX,errQuadW,errQuadRefX,errQuadRefW = \
      DiagUtils.getQuadraturePhysPointsAndWeights(mesh,FemSpace,errorQuadrature,
                                        verbose)
    nquad = Numeric.shape(errQuadX)[1]
    uvals = DiagUtils.getFEMvals(u,errQuadRefX,verbose)
    if verbose > 3:
        print 'errQuadPoints=\n',errQuadX
        print 'numerical solution at errQuadPoints=\n',uvals
    #end if
    eL2  = P1ncMixed.L2errorSFEMvsAF(exactSoln,errQuadX,errQuadW,
                                     uvals,T=0.0)
    print 'L2 error = ',eL2

    outmesh = """mesh%dd """ % nd
    if nd == 1:
        FemSpace.writeFunctionGnuplot(u,'solution')
    elif nd == 2:
        FemSpace.writeFunctionMatlab(u,'solution',outmesh)
    #end nd=2
    #go ahead and plot averages for matlab too? should work as before as long
    #as solution vector is same size as number of triangles
    mesh.writeMeshEnsight(filename='solen')
    FemSpace.writeFunctionEnsight(u,filename='solen',append=False)
    
    # # try to look at solution # #
    x = None
    y = None
    z = None
    aSol = None
    nSol = None

    if nd == 1:
        x = Numeric.arange(nx)/float(nx-1)
        aSol = Numeric.zeros((nx,),Numeric.Float)
        for i in range(nx):
            aSol[i] = exactSoln.uOfXT([x[i],0.,0.0],0.0)
        #end i
        #take advantage of solution layout, vertex dofs are first
        nSol = Numeric.zeros((nx,),Numeric.Float)
        if useCG:
            nSol = u.dof[0:nx]
        else:
            #offset = FemSpace.referenceFiniteElement.localFunctionSpace.dim
            offset = FemSpace.referenceFiniteElement.referenceElement.dim+1
            for i in range(nx-1):
                nSol[i] = u.dof[i*offset]
            #end i
            nSol[-1] = u.dof[-1]
        #end dg
    elif nd == 2:
        x = Numeric.arange(nx)/float(nx-1)
        y = Numeric.arange(ny)/float(ny-1)
        
        #take advantage of solution layout, vertex dofs are first
        if useCG:
            nSol = Numeric.reshape(system.u.dof[0:nx*ny],(nx,ny))
        else:
            print 'fix nSol for 2d DG'
        #end if
        aSol = Numeric.zeros((nx,ny),Numeric.Float)
        for i in range(nx):
            for j in range(ny):
                aSol[i,j] = exactSoln.uOfXT([x[i],y[j],0.0],0.0)
            #end j
        #end i
    elif nd == 3:
        x = Numeric.arange(nx)/float(nx-1)
        y = Numeric.arange(ny)/float(ny-1)
        z = Numeric.arange(nz)/float(nz-1)
        
        #take advantage of solution layout, vertex dofs are first
        nSol = Numeric.reshape(system.u.dof[0:nx*ny*nz],(nx,ny,nz))
        aSol = Numeric.zeros((nx,ny,nz),Numeric.Float)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    aSol[i,j,k] = exactSoln.uOfXT([x[i],y[j],z[k]],0.0)
                #end k
            #end j
        #end i
    #end first if on dim

    import Gnuplot

    if useCG and nd < 3:
        solPlot = Gnuplot.Gnuplot()
        solPlot("set terminal x11")
        aPlot = Gnuplot.Gnuplot()
        aPlot("set terminal x11")

        if nd == 1:

            aPlot.plot(Gnuplot.Data(x,
                                    nSol,
                                    with='linespoints',
                                    title='numerical solution'),
                       Gnuplot.Data(x,
                                    aSol,
                                    with='linespoints',
                                    title='analytical solution'))
            
        elif nd == 2:
        
            solPlot('set parametric')
            solPlot('set data style lines')
            solPlot('set hidden')
            solPlot('set contour base')
            solPlot.xlabel('x')
            solPlot.ylabel('y')
            solPlot.splot(Gnuplot.GridData(nSol,x,y,binary=0,inline=0,
                                            title='numerical solution'))
        
            aPlot('set parametric')
            aPlot('set data style lines')
            aPlot('set hidden')
            aPlot('set contour base')
            aPlot.xlabel('x')
            aPlot.ylabel('y')
            aPlot.splot(Gnuplot.GridData(aSol,x,y,binary=0,inline=0,title='analytical solution'))

            if verbose > 5:            
                print 'u.dofMap.l2g=\n',u.femSpace.dofMap.l2g
            #end verbose

        #end 2d

    #end useCG


    if verbose > 3:
        print 'testLaplacian '
        if nd >= 1:
        #mwf debug
            print 'x=\n',x
        #end if
        if nd >= 2:
            print 'y=\n',y
        #end if
        if nd >= 3:
            print 'z=\n',z
        #end if
        print 'u=\n',u.dof
        print 'exact=\n',aSol
        print 'numerical=\n',nSol
        

    #end verbose > 3
    raw_input('Please press return to continue... \n')



# # # # # # # # # # # # # # # # # # # # # # # # # 
#specific test problems
# # # # # # # # # # # # # # # # # # # # # # # # # 
def collectPoissonDescription(uex,duex,qex,A,A0,f,dir):
    """
    for convenience pull together a Poisson problem description
    consisting of

    uex : exact solution
    duex: derivative of exact solution
    qex : exact flux
    A   : 2 tensor in R^d
    A0  : constant term in tensor description
    f   : source term
    dir : dirichlet bc description
    """
    poisson = {}
    poisson['uex'] = uex
    poisson['duex']= duex
    poisson['qex'] = qex
    poisson['A']   = A
    poisson['A0']  = A0
    poisson['f']   = f
    poisson['DirichletBCs']=dir

    return poisson
#end collect
def uexP0(p,t):
    """
    u = \sum^d_i=1 x_i
    
    """
    #val = p[X] + max(0.,nd-1)*p[Y] + max(0.,nd-2)*p[Z]
    val = p[X]  + p[Y]
    return val
def duexP0(p,t):
    """
    u' = [1,1]
    
    """
    #val = p[X] + max(0.,nd-1)*p[Y] + max(0.,nd-2)*p[Z]
    val = Numeric.ones(2,Numeric.Float)
    return val
def AP0(p):
    """
    conductivity
    """
    return DiagUtils.Ident2
#end def
def qexP0(p,t):
    """
    q = -A\grad u
      = -A*[1,1]^T
    """
    #val = p[X] + max(0.,nd-1)*p[Y] + max(0.,nd-2)*p[Z]
    val = -Numeric.dot(AP0(p),duexP0(p,t))
    return val
def fP0(p):
    """
    f = 0.0
    """
    return 0.0

def dirBCP0(p):
    if p[X] == 0. or p[X] == 1. or p[Y] == 0.0 or p[Y] == 1.0:
        return uexP0
# # # # # # # # # # # # # # # # # # # # # # # # # 
def uexP1(p,t):
    """
    u = \sum^d_i=1 x_i^2
    
    """
    val = p[X]*p[X] + p[Y]*p[Y]
    return val
def duexP1(p,t):
    """
    u' = 2p
    
    """
    val = 2.0*p[0:2]
    return val
def AP1(p):
    """
    conductivity
    """
    return DiagUtils.Ident2
def qexP1(p,t):
    """
    q = -A\grad u
      = -A*[1,1]^T
    """
    #val = p[X] + max(0.,nd-1)*p[Y] + max(0.,nd-2)*p[Z]
    val = -Numeric.dot(AP1(p),duexP1(p,t))
    return val

def fP1(p):
    """
    f = -2*d
    """
    return -4.
def dirBCP1(p):
    if p[X] == 0. or p[X] == 1. or p[Y] == 0.0 or p[Y] == 1.0:
        return uexP1

# # # # # # # # # # # # # # # # # # # # # # # # # 
def uexP3(p,t):
    """
    u = \sum^d_i=1 x_i
    
    """
    #val = p[X] + max(0.,nd-1)*p[Y] + max(0.,nd-2)*p[Z]
    val = p[X]  + p[Y] + p[Z]
    return val
def duexP3(p,t):
    """
    u' = [1,1,1]
    
    """
    val = Numeric.ones(3,Numeric.Float)
    return val
def AP3(p):
    """
    conductivity
    """
    return DiagUtils.Ident3
#end def
def qexP3(p,t):
    """
    q = -A\grad u
      = -A*[1,1]^T
    """
    val = -Numeric.dot(AP3(p),duexP3(p,t))
    return val
def fP3(p):
    """
    f = 0.0
    """
    return 0.0

def dirBCP3(p):
    if p[X] == 0. or p[X] == 1. or p[Y] == 0.0 or p[Y] == 1.0 or p[Z] == 0.0 or p[Z] == 1.0:
        return uexP3
# # # # # # # # # # # # # # # # # # # # # # # # # 
def uexP4(p,t):
    """
    u = \sum^d_i=1 x_i^2
    
    """
    #val = p[X] + max(0.,nd-1)*p[Y] + max(0.,nd-2)*p[Z]
    val = p[X]*p[X]  + p[Y]*p[Y] + p[Z]*p[Z]
    return val
def duexP4(p,t):
    """
    u' = [1,1,1]
    
    """
    val = 2.*p[0:3]
    return val
def AP4(p):
    """
    conductivity
    """
    return DiagUtils.Ident3
#end def
def qexP4(p,t):
    """
    q = -A\grad u
      = -A*[1,1]^T
    """
    val = -Numeric.dot(AP4(p),duexP4(p,t))
    return val
def fP4(p):
    """
    f = -2*d
    """
    return -6.

def dirBCP4(p):
    if p[X] == 0. or p[X] == 1. or p[Y] == 0.0 or p[Y] == 1.0 or p[Z] == 0.0 or p[Z] == 1.0:
        return uexP4

# # # # # # # # # # # # # # # # # # # # # # # # # 
def uexP5(p,t):
    """
    u = \sum^d_i=1 x_i^2
    d = 1
    """
    val = p[X]*p[X]
    return val
def duexP5(p,t):
    """
    u' = 2p
    
    """
    val = 2.0*p[X]
    return val
def AP5(p):
    """
    conductivity
    """
    return DiagUtils.Ident1
def qexP5(p,t):
    """
    q = -A\grad u
      = -A*[1,1]^T
    """
    #val = p[X] + max(0.,nd-1)*p[Y] + max(0.,nd-2)*p[Z]
    val = -Numeric.dot(AP5(p),duexP5(p,t))
    return val

def fP5(p):
    """
    f = -2*d
    """
    return -2.
def dirBCP5(p):
    if p[X] == 0. or p[X] == 1.: 
        return uexP5
# # # # # # # # # # # # # # # # # # # # # # # # # 

## @}
