from math import *
import cPickle
from MeshTools import *
from FemTools import *
from QuadTools import *
from LinearSolvers import *
from NonlinearSolvers import *
import Gnuplot
from NormTools import *
from AnalyticalSolutions import *
from HamiltonJacobi import *
"""
A module for testing Hamilton-Jacobi models.
"""

## \defgroup HamiltonJacobiTests HamiltonJacobiTests
#
# A module for testing Hamilton-Jacobi models.  
#
# @{

def runTests():
    #First define the mesh independent problem definition
    #consisting of the following information
    #see examples below for definitions
    nd={}  #number of spatial dimensions
    getDirichletConditions={}
    coefficients={}
    analyticalSolution={}
    timeIntegration = {}
    getInitialConditions ={}
    T = {}
    fullNewtonFlag = {}
    headerFile = open('reportHeader.tex','r')
    header = headerFile.read()
    report = open('HamiltonJacobiReport.tex','w')
    report.write(header)
    #Define test problems
    testProblems = []
    #1D
    #
    #  h  = (1,0,0) * grad(u)
    #
    #  u(x,0) is the signed distance function
    #  from circle of radius 1/8 centered at 1/4,1/2
    #
    test='ConstantVelocityCircle'
    testProblems.append(test)
    timeIntegration[test]=ForwardEuler
    class ConstantVelocityLevelSet(HamiltonJacobiCoefficients):
        def __init__(self,B):
            self.B = B
        def evaluate(self,
                     t,
                     x,
                     u,grad_u,
                     m,dm,
                     h,dh,
                     rh):
            rh[:]=0.0
            for i in range(len(h)):
                h[i]=Numeric.dot(self.B,grad_u[i])
                dh[i][:]=self.B
            m[:] = u
            dm[:] = 1.0
    nd[test]=2
    N=3.0
    b = Numeric.array([0.0,-1.0])
    coefficients[test]=ConstantVelocityLevelSet(b)
    class ConstantVelocityCircle:
        def __init__(self,radius,B,startX=0.25,startY=0.5):
            self.radius = radius
            self.B = B
            self.startX=startX
            self.startY=startY
        def uOfXT(self,x,t):
            centerX = self.B[0]*t + self.startX 
            centerY = self.B[1]*t + self.startY
            return sqrt((x[0]-centerX)**2 + (x[1] - centerY)**2) - self.radius
    T[test]=0.5
#      cvelSol = ConstantVelocityCircle(1.0/8.0,b)
    cvelSol = ConstantVelocityCircle(1.0/8.0,b,startX=0.5,startY=0.5)
    analyticalSolution[test] = cvelSol
    def get_cvelSolDBC(x):
        pass
#          if (x[X] == 0.0 or
#              x[X] == 1.0 or
#              x[Y] == 0.0 or
#              x[Y] == 1.0):
#              return cvelSol.uOfXT
    getDirichletConditions[test]=get_cvelSolDBC
    getInitialConditions[test] = analyticalSolution[test]
    fullNewtonFlag[test] = False#True
    coefficients[test].mass = 'linear'
    coefficients[test].hamiltonian = 'linear'
    #
    #  h  = (1,0,0) * grad(u)
    #
    #  u(x,0) is the signed distance function
    #  from circle of radius 1/8 centered at 1/4,1/2
    #
    test='RotatingCircle'
    testProblems.append(test)
    timeIntegration[test]=ForwardEuler
    class Rotating2DVelocityLevelSet(HamiltonJacobiCoefficients):
        def __init__(self):
            pass
        def evaluate(self,
                     t,
                     x,
                     u,grad_u,
                     m,dm,
                     h,dh,
                     rh):
            rh[:]=0.0
            for i in range(len(h)):
                vx = 2*pi*(x[i][1] - 0.5)
                vy = 2*pi*(0.5     - x[i][0]) 
                h[i]=(vx*grad_u[i][0] +
                      vy*grad_u[i][1])
                dh[i][0]=vx
                dh[i][1]=vy
            m[:] = u
            dm[:] = 1.0
    nd[test]=2
    N=3.0
    coefficients[test]=Rotating2DVelocityLevelSet()
    class Rotating2DVelocityCircle:
        def __init__(self,radius):
            self.radius = radius
        def uOfXT(self,x,t):
            centerX = 0.25*sin(2*pi*t) + 0.5
            centerY = 0.25*cos(2*pi*t) + 0.5
            return sqrt((x[0]-centerX)**2 + (x[1] - centerY)**2) - self.radius
    T[test]=0.5
    analyticalSolution[test] = Rotating2DVelocityCircle(1.0/8.0)
    def dbc(n,t):
        return analyticalSolution[test].uOfXT(n.p,t)
    def getSolutionDirichletConditions(n):
        pass
#              if (n.p[X] == 0.0 or
#                  n.p[X] == 1.0 or
#                  n.p[Y] == 0.0 or
#                  n.p[Y] == 1.0):
#                  return dbc
#              pass
    getDirichletConditions[test]=getSolutionDirichletConditions
    getInitialConditions[test] = analyticalSolution[test]
    fullNewtonFlag[test] = False#True
    coefficients[test].mass = 'linear'
    coefficients[test].hamiltonian = 'linear'
    #
    #  h  = |\grad u\|
    #
    #  u(x,0) is the signed distance function
    #  from the circle of radius 1/8 centered at 1/2,1/2
    #
    test='ExpandingCircle'
    testProblems.append(test)
    timeIntegration[test]=ForwardEuler
    class UnitNormalVelocity(HamiltonJacobiCoefficients):
        def __init__(self):
            pass
        def evaluate(self,
                     t,
                     x,
                     u,grad_u,
                     m,dm,
                     h,dh,
                     rh):
            rh[:]=0.0
            for i in range(len(h)):
                h[i]=sqrt(grad_u[i][0]**2 +
                          grad_u[i][1]**2)
                if h[i] != 0.0:
                    dh[i][0]=grad_u[i][0]/h[i]
                    dh[i][1]=grad_u[i][1]/h[i]
                else:
                    dh[i][0] = grad_u[i][0]/1.0e-8
                    dh[i][1] = grad_u[i][1]/1.0e-8
            m[:] = u
            dm[:] = 1.0
    nd[test]=2
    N=3.0
    coefficients[test]=UnitNormalVelocity()
    class ExpandingCircle:
        def __init__(self,radius):
            self.radius = radius
        def uOfXT(self,x,t):
            centerX = 0.5
            centerY = 0.5
            return sqrt((x[0]-centerX)**2 + (x[1] - centerY)**2) - (self.radius + t) 
    analyticalSolution[test] = ExpandingCircle(1.0/8.0)
    T[test]=1.0-(1.0/2.0 + 1.0/8.0)
    def dbc(n,t):
        return analyticalSolution[test].uOfXT(n.p,t)
    def getSolutionDirichletConditions(n):
        pass
#              if (n.p[X] == 0.0 or
#                  n.p[X] == 1.0 or
#                  n.p[Y] == 0.0 or
#                  n.p[Y] == 1.0):
#                  return dbc
#              pass
    getDirichletConditions[test]=getSolutionDirichletConditions
    getInitialConditions[test] = analyticalSolution[test]
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'linear'
    coefficients[test].hamiltonian = 'nonlinear'
    #
    #  h  = |\grad u\|
    #
    #  u(x,0) is the signed distance function
    #  from the sqaure with side 2*1/8 centered at 1/2,1/2
    #
    test='ExpandingSquare'
    testProblems.append(test)
    timeIntegration[test]=ForwardEuler
    nd[test]=2
    N=3.0
    coefficients[test]=UnitNormalVelocity()
    class ExpandingSquare:
        def __init__(self,radius):
            self.radius = radius
        def uOfXT(self,x,t):
            centerX = 0.5
            centerY = 0.5
            oX = x[0] - centerX
            oY = x[1] - centerY
            if oY < oX:
                if oY > -oX:
                    return oX - (self.radius + t)
                else:
                    return -oY - (self.radius + t)
            else:
                if oY > -oX:
                    return oY - (self.radius + t)
                else:
                    return -oX - (self.radius + t)
    analyticalSolution[test] = ExpandingSquare(1.0/8.0)
    T[test]=1.0-(1.0/2.0 + 1.0/8.0)
    def dbc(n,t):
        return analyticalSolution[test].uOfXT(n.p,t)
    def getSolutionDirichletConditions(n):
        pass
#              if (n.p[X] == 0.0 or
#                  n.p[X] == 1.0 or
#                  n.p[Y] == 0.0 or
#                  n.p[Y] == 1.0):
#                  return dbc
    getDirichletConditions[test]=getSolutionDirichletConditions
    getInitialConditions[test] = analyticalSolution[test]
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'linear'
    coefficients[test].hamiltonian = 'nonlinear'
    #
    #  h  = |\grad u\|
    #
    #  u(x,0) is the signed distance function
    #  from the sqaure with side 2*1/8 centered at 1/2,1/2
    #
    test='MergingCircles'
    testProblems.append(test)
    timeIntegration[test]=ForwardEuler
    nd[test]=2
    N=3.0
    coefficients[test]=UnitNormalVelocity()
    class MergingCircles:
        def __init__(self,radius):
            self.radius = radius
        def uOfXT(self,x,t):
            centerX1 = 0.5 - 1.0/8.0 - 1.0/16.0
            centerX2 = 0.5 + 1.0/8.0 + 1.0/16.0
            centerY = 0.5
            oX1 = x[0] - centerX1
            oX2 = x[0] - centerX2
            oY = x[1] - centerY
            r1 = sqrt((x[0]-centerX1)**2 + (x[1] - centerY)**2)
            r2 = sqrt((x[0]-centerX2)**2 + (x[1] - centerY)**2)
            if r1 < r2:
                return r1 - (self.radius + t)
            else:
                return r2 - (self.radius + t)
    analyticalSolution[test] = MergingCircles(1.0/8.0)
    T[test]=1.0-(1.0/2.0 + 1.0/8.0)
    def dbc(n,t):
        return analyticalSolution[test].uOfXT(n.p,t)
    def getSolutionDirichletConditions(n):
        pass
#              if (n.p[X] == 0.0 or
#                  n.p[X] == 1.0 or
#                  n.p[Y] == 0.0 or
#                  n.p[Y] == 1.0):
#                  return dbc
    getDirichletConditions[test]=getSolutionDirichletConditions
    getInitialConditions[test] = analyticalSolution[test]
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'linear'
    coefficients[test].hamiltonian = 'nonlinear'
    #
    #  h  = |\grad u\| - 1
    #
    #  u(x,0) is the signed distance function squared
    #  from the circle of radius 1/8 centered at 1/2,1/2
    #
    test='MergineCircleEikonal'
    testProblems.append(test)
    timeIntegration[test]=ForwardEuler#NoIntegration
    class UnitNormalVelocityMergingCirclesEikonal(HamiltonJacobiCoefficients):
        def __init__(self):
            self.mass='linear'
            self.hamiltonian='nonlinear'
        def initializeSignFunction(self,uq,grad_uq,diameter):
            self.Sq=Numeric.zeros(uq.shape,
                                  Numeric.Float)
        def updateSignFunction(self,uq,grad_uq,diameter):
            for i in range(uq.shape[0]):
                norm_grad_u_2=0.0
                for I in range(grad_uq.shape[1]):
                    norm_grad_u_2+=grad_uq[i,I]*grad_uq[i,I]
#                 self.Sq[i] =  uq[i]/(math.sqrt(uq[i]*uq[i] + norm_grad_u_2*4.0*diameter*diameter))
                self.Sq[i] =  uq[i]/(math.sqrt(uq[i]*uq[i] + diameter*diameter))
                #print self.Sq[i]
#                  if  uq[i] > 0.0:
#                      self.Sq[i]=1.0
#                  elif uq[i] < 0.0:
#                      self.Sq[i]=-1.0
#                  else:
#                      self.Sq[i]=0.0
        def evaluate(self,
                     t,
                     x,
                     u,grad_u,
                     m,dm,
                     h,dh,
                     rh):
            rh[:]=-1.0
            dm[:]=1.0
            h[:]=0.0
            for i in range(u.shape[0]):
                m[i]=u[i]
                for I in range(grad_u.shape[1]):
                    h[i]+=grad_u[i,I]*grad_u[i,I]
                h[i]=math.sqrt(h[i])
                for  I in range(grad_u.shape[1]):
                    if  h[i] != 0.0:
                        dh[i,I]=(self.Sq[i]*grad_u[i,I])/h[i]
                    else:
                        dh[i,I]=(self.Sq[i]*grad_u[i,I])/1.0e-8
                h[i]+=rh[i]
                h[i]*=self.Sq[i]
                rh[i]*=self.Sq[i]
#          def __init__(self):
#              pass
#          def evaluate(self,
#                       t,
#                       x,
#                       u,grad_u,
#                       m,dm,
#                       h,dh,
#                       rh):
#              rh[:]=-1.0
#              dm[:]=1.0
#              m[:]=u
#              h[:]=0.0
#              for i in range(u.shape[0]):
#                  for I in range(grad_u.shape[1]):
#                      h[i]+=grad_u[i,I]*grad_u[i,I]
#                  h[i]=math.sqrt(h[i])
#                  for  I in range(grad_u.shape[1]):
#                      if  h[i] != 0.0:
#                          dh[i,I]=grad_u[i,I]/h[i]
#                      else:
#                          dh[i,I]=grad_u[i,I]/1.0e-8
#                  h[i]+=rh[i]
    nd[test]=2
    N=3.0
    coefficients[test]=UnitNormalVelocityMergingCirclesEikonal()
    class MergingCirclesEikonal:
        def __init__(self,radius):
            self.radius = radius
        def uOfXT(self,x,t):
            centerX1 = 0.5 -2.0*self.radius - 1.0/8.0 - 1.0/16.0
            centerX2 = 0.5 #+ 1.0/8.0 + 1.0/16.0
#                  centerX1 = 0.5 
#                  centerX2 = 0.5 
            centerY = 0.5
            oX1 = x[0] - centerX1
            oX2 = x[0] - centerX2
            oY = x[1] - centerY
            r1 = sqrt((x[0]-centerX1)**2 + (x[1] - centerY)**2)
            r2 = sqrt((x[0]-centerX2)**2 + (x[1] - centerY)**2)
            if r1 < r2:
                return (r1 - self.radius)**3#*abs(r1 - self.radius)
            else:
                return (r2 - self.radius)**3#*abs(r2 - self.radius)
#             return (math.sqrt(x[1]**2 + x[0]**2)-0.75)**3
#             return -(x[1] - (x[0]-0.5)**2 -0.5 )**3
    mceSol = MergingCirclesEikonal(1.0/8.0)
    analyticalSolution[test] = mceSol
    T[test]=10.0
    def getSolutionDirichletConditions(x):
            pass
    getDirichletConditions[test]=getSolutionDirichletConditions
    getInitialConditions[test] = mceSol
    fullNewtonFlag[test] = True
    coefficients[test].mass = 'linear'
    coefficients[test].hamiltonian = 'nonlinear'
    for test in testProblems[:1]:
        computeEigenvalues=False
        if computeEigenvalues:
            linearSolverType= levelLinearSolverType = 'DenseLU'
            levelNonlinearSolverType = 'Newton'
            nonlinearSolverType = 'NLNI'
        else:
            linearSolverType = levelLinearSolverType = 'SparseLU'
            #linearSolverType= levelLinearSolverType = 'DenseLU' 
            #linearSolverType= levelLinearSolverType = 'Jacobi' 
            #linearSolverType= levelLinearSolverType = 'GaussSeidel' 
            #linearSolverType= levelLinearSolverType = 'StarILU' 
            #linearSolverType = 'NI'
            #levelLinearSolverType = 'MGM'
            #smootherType = 'StarILU'
            #smootherType = 'GaussSeidel'
            #smootherType = 'Jacobi'
            levelNonlinearSolverType = 'Newton'
            #levelNonlinearSolverType = 'NLJacobi'
            #levelNonlinearSolverType = 'NLGaussSeidel'
            #levelNonlinearSolverType = 'NLStarILU'
            #levelNonlinearSolverType = 'FAS'
            nonlinearSolverType = 'NLNI'
            #nonlinearSolverType = 'Newton'
        tolFac = 0.01
        linTolFac = 0.001
        runCFL = 0.1
        DG=False
        stabilization='2'
        shockCapturing=None#'2'
        shockCapturingDiffusion=0.5
        quadratureOrder=3
        preSmooths = 2
        postSmooths = 2
        cycles = 2
        nLevels=5
        computeSpaceTimeError=False
        nn=3 #number of nodes on the coarsest mesh
        print "Starting Test "+`test`
        print "Setting up quadrature"
        quadrature={}
        gq = SimplexGaussQuadrature(nd[test])
        #gq = SimplexLobattoQuadrature(nd[test])
        gq.setOrder(quadratureOrder)
        for integral in OneLevelHamiltonJacobi.integralKeys:
            quadrature[integral] = gq 
        #
        #define the mesh hierarchy
        #
        print "Setting up MultilevelMesh"
        mlMesh = None
        mlMeshFileName = "mlMesh%dD.%d" % (nd[test],nLevels)
        try:
            mlMeshFile = open(mlMeshFileName,'rb')
            print "reading mesh"
            mlMesh = cPickle.load(mlMeshFile)
            print "done reading mesh"
        except:
            print "generating mesh"
            mlMeshFile = open(mlMeshFileName,'wb')
            if nd[test]==1:
                mlMesh = MultilevelEdgeMesh(nn,1,1,refinementLevels=nLevels)
            elif nd[test]==2:
                mlMesh = MultilevelTriangularMesh(nn,nn,1,
                                                  refinementLevels=nLevels)
            elif nd[test]==3:
                mlMesh = MultilevelTetrahedralMesh(nn,nn,nn,
                                                   refinementLevels=nLevels)
            cPickle.dump(mlMesh,mlMeshFile,protocol=cPickle.HIGHEST_PROTOCOL)
            print "done generating mesh"
        print "Setting up MultilevelHamiltonJacobi"
        tolList=[]
        linTolList=[]
        for l in range(nLevels):
            mlMesh.meshList[l].computeGeometricInfo()
            tolList.append(tolFac*(mlMesh.meshList[l].h**2))
            linTolList.append(linTolFac*(mlMesh.meshList[l].h**2))
        atol = min(tolList)
        lin_atol = min(linTolList)
        if computeEigenvalues or linearSolverType == 'DenseLU':
            MatType = Mat
        else:
            MatType = SparseMat
        mlHamiltonJacobi = MultilevelHamiltonJacobi(
            nd[test],
            mlMesh,
            C0_AffineLinearOnSimplexWithNodalBasis,
            C0_AffineLinearOnSimplexWithNodalBasis,
            MatType,
            getDirichletConditions[test],
            coefficients[test],
            quadrature,
            stabilization,
            shockCapturing,
            shockCapturingDiffusion,
            timeIntegration[test])
        for l in range(nLevels):
            mlHamiltonJacobi.modelList[l].coefficients.initializeSignFunction(mlHamiltonJacobi.modelList[l].q['u'].flat,
                                                                               Numeric.reshape(mlHamiltonJacobi.modelList[l].q['grad(u)'].flat,
                                                                                               (mlHamiltonJacobi.modelList[l].nQuadraturePoints_global,
                                                                                                mlHamiltonJacobi.modelList[l].nSpace_global)),
                                                                              1.5*mlHamiltonJacobi.modelList[l].mesh.h)
        utmp = []
        uqtmp = Numeric.zeros(mlHamiltonJacobi.modelList[-1].q['u'].shape,Numeric.Float)
        asolq = Numeric.zeros(mlHamiltonJacobi.modelList[-1].q['u'].shape,Numeric.Float)
        for l in range(nLevels):
            utmp.append(FiniteElementFunction(mlHamiltonJacobi.modelList[l].trialSpace))
        print "Setting up LinearSolver"
        levelLinearSolver = None
        if levelLinearSolverType == 'MGM':
            prolongList=[]
            restrictList=[]
            preSmootherList=[]
            postSmootherList=[]
            mgItsList=[]
            prolongList = mlHamiltonJacobi.meshTransfers.prolongList
            restrictList = mlHamiltonJacobi.meshTransfers.restrictList
            for l in range(nLevels):
                mgItsList.append(cycles)
                if l > 0:
                    if smootherType == 'Jacobi':
                        preSmootherList.append(
                            Jacobi(L=mlHamiltonJacobi.jacobianList[l],
                                   weight=4.0/5.0,
                                   maxIts=preSmooths,
                                   convergenceTest = 'its',
                                   computeRates = True,
                                   printInfo = False))
                        postSmootherList.append(
                            Jacobi(L=mlHamiltonJacobi.jacobianList[l],
                                   weight=4.0/5.0,
                                   maxIts=postSmooths,
                                   convergenceTest = 'its',
                                   computeRates = True,
                                   printInfo = False))
                    elif smootherType == 'GaussSeidel':
                        preSmootherList.append(
                            GaussSeidel(connectionList = \
                                        mlHamiltonJacobi.modelList[l].freeNodeStarList,
                                        L=mlHamiltonJacobi.jacobianList[l],
                                        weight=1.0,
                                        maxIts =  preSmooths,
                                        convergenceTest = 'its',
                                        computeRates = True,
                                        printInfo = False))
                        postSmootherList.append(
                            GaussSeidel(connectionList = \
                                        mlHamiltonJacobi.modelList[l].freeNodeStarList,
                                        L=mlHamiltonJacobi.jacobianList[l],
                                        weight=1.0,
                                        maxIts =  postSmooths,
                                        convergenceTest = 'its',
                                        computeRates = True,
                                        printInfo = False))
                    elif smootherType == 'StarILU':
                        preSmootherList.append(
                            StarILU(connectionList = \
                                    mlHamiltonJacobi.modelList[l].freeNodeStarList,
                                    L=mlHamiltonJacobi.jacobianList[l],
                                    weight=1.0,
                                    maxIts =  preSmooths,
                                    convergenceTest = 'its',
                                    computeRates = True,
                                    printInfo = False))
                        postSmootherList.append(
                            StarILU(connectionList = \
                                    mlHamiltonJacobi.modelList[l].freeNodeStarList,
                                    L=mlHamiltonJacobi.jacobianList[l],
                                    weight=1.0,
                                    maxIts =  postSmooths,
                                    convergenceTest = 'its',
                                    computeRates = True,
                                    printInfo = False))
                    else:
                        print "smootherType unrecognized"
                else:
                    preSmootherList.append([])
                    postSmootherList.append([])
                    coarseSolver = SparseLU(
                        L=mlHamiltonJacobi.jacobianList[l])
            levelLinearSolver = MGM(prolongList = prolongList,
                                    restrictList = restrictList,
                                    LList = mlHamiltonJacobi.jacobianList,
                                    preSmootherList = preSmootherList,
                                    postSmootherList = postSmootherList,
                                    coarseSolver = coarseSolver,
                                    mgItsList = mgItsList)
            levelLinearSolverList = levelLinearSolver.solverList
        elif levelLinearSolverType == 'DenseLU':
            levelLinearSolverList = []
            for l in range(nLevels):
                levelLinearSolverList.append(
                    DenseLU(mlHamiltonJacobi.jacobianList[l]))
            levelLinearSolver = levelLinearSolverList
        elif linearSolverType == 'SparseLU':
            levelLinearSolverList=[]
            for l in range(nLevels):
                levelLinearSolverList.append(
                    SparseLU(mlHamiltonJacobi.jacobianList[l]))
            levelLinearSolver = levelLinearSolverList
        elif linearSolverType == 'Jacobi':
            levelLinearSolverList=[]
            for l in range(nLevels):
                levelLinearSolverList.append(
                    Jacobi(L=mlHamiltonJacobi.jacobianList[l],
                           weight=4.0/5.0,
                           maxIts=500,
                           convergenceTest = 'r',
                           rtol_r = linTolList[l],
                           atol_r = lin_atol,
                           computeRates = True,
                           printInfo = False))
            levelLinearSolver = levelLinearSolverList
        elif linearSolverType == 'GaussSeidel':
            levelLinearSolverList=[]
            for l in range(nLevels):
                levelLinearSolverList.append(
                    GaussSeidel(connectionList = \
                                mlHamiltonJacobi.modelList[l].freeNodeStarList,
                                L=mlHamiltonJacobi.jacobianList[l],
                                weight=1,
                                maxIts=500,
                                convergenceTest = 'r',
                                rtol_r = linTolList[l],
                                atol_r = lin_atol,
                                computeRates = True,
                                printInfo = False))
            levelLinearSolver = levelLinearSolverList
        elif linearSolverType == 'StarILU':
            levelLinearSolverList=[]
            for l in range(nLevels):
                levelLinearSolverList.append(
                    StarILU(connectionList = \
                            mlHamiltonJacobi.modelList[l].freeNodeStarList,
                            L=mlHamiltonJacobi.jacobianList[l],
                            weight=1.0,
                            maxIts=500,
                            convergenceTest = 'r',
                            rtol_r = linTolList[l],
                            atol_r = lin_atol,
                            computeRates = True,
                            printInfo = False))
            levelLinearSolver = levelLinearSolverList
        else:
            print "Unknown level linear solver "+ levelLinearSolverType
        linearSolver = None
        if linearSolverType == 'NI':
            ni = NI(solverList = levelLinearSolverList,
                    prolongList = prolongList,
                    restrictList = restrictList,
                    maxIts  = 500,
                    tolList = linTolList,
                    atol    = lin_atol,
                    printInfo=False)
            linearSolver = []
            #nested iteration knows how to solve any of the levels
            for l in range(nLevels):
                linearSolver.append(ni)
        elif (linearSolverType == 'DenseLU' or
              linearSolverType == 'SparseLU' or
              linearSolverType == 'Jacobi' or
              linearSolverType == 'GaussSeidel' or
              linearSolverType == 'StarILU'):
            linearSolver = levelLinearSolver
            for l in range(nLevels):
                linearSolver[l].printInfo=False
        else:
            print "Unknown linear solver %s" % linearSolverType
        if (linearSolverType == 'DenseLU' or
            linearSolverType == 'SparseLU'):
            directSolverFlag=True
        else:
            directSolverFlag=False
        print "Setting up NonlinearSolver"
        levelNonlinearSolverList=[]
        if levelNonlinearSolverType == 'Newton':
            for l in range(nLevels):
                levelNonlinearSolverList.append(
                    Newton(linearSolver=linearSolver[l],
                           F=mlHamiltonJacobi.modelList[l],
                           J=mlHamiltonJacobi.jacobianList[l],
                           rtol_r=tolList[l],
                           atol_r=atol,
                           maxIts=500,
                           convergenceTest = 'r',
                           printInfo=True,#False,
                           fullNewton=fullNewtonFlag[test],
                           directSolver=directSolverFlag))
        elif levelNonlinearSolverType == 'NLJacobi':
            for l in range(nLevels):
                levelNonlinearSolverList.append(
                    NLJacobi(F=mlHamiltonJacobi.modelList[l],
                             J=mlHamiltonJacobi.jacobianList[l],
                             rtol_r=tolList[l],
                             atol_r=atol,
                             maxIts=500,
                             convergenceTest = 'r',
                             weight=4.0/5.0,
                             printInfo=False,
                             fullNewton=fullNewtonFlag[test]))
        elif levelNonlinearSolverType == 'NLGaussSeidel':
            for l in range(nLevels):
                levelNonlinearSolverList.append(
                    NLGaussSeidel(F=mlHamiltonJacobi.modelList[l],
                                  J=mlHamiltonJacobi.jacobianList[l],
                                  connectionList = \
                                  mlHamiltonJacobi.modelList[l].freeNodeStarList,
                                  rtol_r=tolList[l],
                                  atol_r=atol,
                                  maxIts=500,
                                  convergenceTest = 'r',
                                  weight=1.0,
                                  printInfo=False,
                                  fullNewton=fullNewtonFlag[test]))
        elif levelNonlinearSolverType == 'NLStarILU':
            for l in range(nLevels):
                levelNonlinearSolverList.append(
                    NLStarILU(F=mlHamiltonJacobi.modelList[l],
                              J=mlHamiltonJacobi.jacobianList[l],
                              connectionList = \
                              mlHamiltonJacobi.modelList[l].freeNodeStarList,
                              rtol_r=tolList[l],
                              atol_r=atol,
                              maxIts=500,
                              convergenceTest = 'r',
                              weight=1.0,
                              printInfo=False,
                              fullNewton=fullNewtonFlag[test]))
        elif levelNonlinearSolverType == 'FAS':
            prolongList=[]
            restrictList=[]
            restrictionRowSumList=[]
            preSmootherList=[]
            postSmootherList=[]
            mgItsList=[]
            prolongList = mlHamiltonJacobi.meshTransfers.prolongList
            restrictList = mlHamiltonJacobi.meshTransfers.restrictList
            restrictionRowSumList = mlHamiltonJacobi.meshTransfers.restrictSumList
            for l in range(nLevels):
                mgItsList.append(cycles)
                if l > 0:
                    if smootherType == 'Jacobi':
                        preSmootherList.append(
                            NLJacobi(F=mlHamiltonJacobi.modelList[l],
                                     J=mlHamiltonJacobi.jacobianList[l],
                                     weight=4.0/5.0,
                                     maxIts=preSmooths,
                                     convergenceTest='its',
                                     printInfo=False,
                                     fullNewton=fullNewtonFlag[test]))
                        postSmootherList.append(
                            NLJacobi(F=mlHamiltonJacobi.modelList[l],
                                     J=mlHamiltonJacobi.jacobianList[l],
                                     weight=4.0/5.0,
                                     maxIts=postSmooths,
                                     convergenceTest='its',
                                     printInfo=False,
                                     fullNewton=fullNewtonFlag[test]))
                    elif smootherType == 'GaussSeidel':
                        preSmootherList.append(
                            NLGaussSeidel(connectionList = \
                                          mlHamiltonJacobi.modelList[l].freeNodeStarList,
                                          F=mlHamiltonJacobi.modelList[l],
                                          J=mlHamiltonJacobi.jacobianList[l],
                                          weight=1.0,
                                          maxIts=preSmooths,
                                          convergenceTest='its',
                                          printInfo=False,
                                          fullNewton=fullNewtonFlag[test]))
                        postSmootherList.append(
                            NLGaussSeidel(connectionList = \
                                          mlHamiltonJacobi.modelList[l].freeNodeStarList,
                                          F=mlHamiltonJacobi.modelList[l],
                                          J=mlHamiltonJacobi.jacobianList[l],
                                          weight=1.0,
                                          maxIts=postSmooths,
                                          convergenceTest='its',
                                          printInfo=False,
                                          fullNewton=fullNewtonFlag[test]))
                    elif smootherType == 'StarILU':
                        preSmootherList.append(
                            NLStarILU(connectionList = \
                                      mlHamiltonJacobi.modelList[l].freeNodeStarList,
                                      F=mlHamiltonJacobi.modelList[l],
                                      J=mlHamiltonJacobi.jacobianList[l],
                                      weight=1.0,
                                      maxIts=preSmooths,
                                      convergenceTest='its',
                                      printInfo=False,
                                      fullNewton=fullNewtonFlag[test]))
                        postSmootherList.append(
                            NLStarILU(connectionList = \
                                      mlHamiltonJacobi.modelList[l].freeNodeStarList,
                                      F=mlHamiltonJacobi.modelList[l],
                                      J=mlHamiltonJacobi.jacobianList[l],
                                      weight=1.0,
                                      maxIts=postSmooths,
                                      convergenceTest='its',
                                      printInfo=False,
                                      fullNewton=fullNewtonFlag[test]))
                    else:
                        print "smootherType unrecognized"
                else:
                    preSmootherList.append([])
                    postSmootherList.append([])
                    coarseSolver = Newton(F=mlHamiltonJacobi.modelList[l],
                                          J= \
                                        mlHamiltonJacobi.jacobianList[l],
                                          linearSolver= \
                                          SparseLU(
                                        mlHamiltonJacobi.jacobianList[l]),
                                          rtol_r=tolList[0],
                                          atol_r=atol,
                                          convergenceTest='r',
                                          maxIts=500,
                                          printInfo=False,
                                          fullNewton=fullNewtonFlag[test],
                                          directSolver=True)
            levelNonlinearSolver = FAS(prolongList = prolongList,
                                       restrictList = restrictList,
                                       restrictSumList = restrictionRowSumList,
                                       FList = mlHamiltonJacobi.modelList,
                                       preSmootherList = preSmootherList,
                                       postSmootherList = postSmootherList,
                                       coarseSolver = coarseSolver,
                                       mgItsList = mgItsList,
                                       printInfo=False)
            levelNonlinearSolverList = levelNonlinearSolver.solverList
        else:
            print "Unknown level nonlinear solver "+ \
                   levelNonlinearSolverType
        if nonlinearSolverType == 'NLNI':
            nonlinearSolver = NLNI(fList = mlHamiltonJacobi.modelList,
                                   solverList = levelNonlinearSolverList,
                                   prolongList = \
                        mlHamiltonJacobi.meshTransfers.prolong_bcList,
                                   restrictList = \
                        mlHamiltonJacobi.meshTransfers.restrict_bcList,
                                   restrictSumList = \
                        mlHamiltonJacobi.meshTransfers.restrict_bcSumList,
                                   maxIts = 500,
                                   tolList = tolList,
                                   atol=atol,
                                   printInfo=True)
        elif nonlinearSolverType == 'Newton':
            nonlinearSolver = levelNonlinearSolverList
        else:
            print "Unknown nonlinearSolverType " + nonlinearSolverType
        print "Running Solver"
        if timeIntegration[test] == NoIntegration:
            T[test] = 1.0
        tn = 0.0
        nSteps = 0
        if getInitialConditions[test] != None:
            mlHamiltonJacobi.setInitialConditions(getInitialConditions[test],tn)
        nx=[]
        ny=[]
        x=[]
        y=[]
        aSol=[]
        solPlot = Gnuplot.Gnuplot()
        solPlot("set terminal x11")
        aSolPlot = Gnuplot.Gnuplot()
        aSolPlot("set terminal x11")
        for l in range(nLevels):
            if DG != True:
                if nd[test]==1:
                    solPlot.title(test)
                    aSolPlot.title(test)
                    nap=101
                    dxap=Numeric.array([1.0/(nap - 1.0),0.0,0.0])
                    P = [(i*dxap) for i in range(nap)]
                    Px = [x[0] for x in P]
                    solPlot.plot(Gnuplot.Data(mlMesh.meshList[l].nodeArray[:,0],
                                              mlHamiltonJacobi.modelList[l].u.dof,
                                              with='linespoints',
                                              title='numerical solution (initial condition/guess)'))
                    if analyticalSolution[test] != None:
                        aSolPlot.plot(Gnuplot.Data(Px,
                                                   [analyticalSolution[test].uOfXT(x,tn) for x in P],
                                                   with='lines',
                                                   title='analytical solution'))
                    else:
                        aSolPlot=solPlot
                elif nd[test]==2:
                    nx.append((nn-1)*(2**l)+1)
                    ny.append(nx[l])
                    x.append(Numeric.arange(nx[l])/float(nx[l]-1))
                    y.append(Numeric.arange(nx[l])/float(nx[l]-1))
                    nSol = Numeric.reshape(mlHamiltonJacobi.modelList[l].u.dof,(nx[l],ny[l]))
                    solPlot('set parametric')
                    solPlot('set data style lines')
                    solPlot('set hidden')
                    solPlot('set contour base')
                    #solPlot('set cntrparam levels incremental 0.0,0.1,1.0')
                    solPlot.xlabel('x')
                    solPlot.ylabel('y')
                    solPlot.splot(Gnuplot.GridData(nSol,
                                                   x[l],
                                                   y[l],
                                                   binary=0,
                                                   inline=0
                                                   ))
                    if analyticalSolution[test] != None:
                        aSol.append(Numeric.zeros((nx[l],ny[l]),Numeric.Float))
                        for i in range(nx[l]):
                            for j in range(ny[l]):
                                aSol[l][i,j] = analyticalSolution[test].uOfXT(Numeric.array([x[l][i],y[l][j],0.0]),tn)
                        aSolPlot('set parametric')
                        aSolPlot('set data style lines')
                        aSolPlot('set hidden')
                        aSolPlot('set contour base')
                        #aSolPlot('set cntrparam levels incremental 0.0,0.1,1.0')
                        aSolPlot.xlabel('x')
                        aSolPlot.ylabel('y')
                        aSolPlot.splot(Gnuplot.GridData(aSol[l],
                                                       x[l],
                                                       y[l],
                                                       binary=0,
                                                       inline=0
                                       ))
        mlHamiltonJacobi.modelList[-1].timeIntegration.runCFL = runCFL
        if int(T[test]/mlMesh.meshList[-1].h**2) != 0.0:
            DTSET = 0.01*T[test]/int(T[test]/mlMesh.meshList[-1].h**2)
            #DTSET = 0.01*T[test]/int(T[test]/mlMesh.meshList[-1].h)
        else:
            DTSET = T[test]
        DTSET=None
        eSpace = {}
        eSpaceTime = {}
        import sys
        tstring=None
        eraseTime='\b\b\b\b\b\b\b\b\b\b\b\b'
        for l in range(nLevels):
            mlHamiltonJacobi.modelList[l].updateZeroLevelSetDirichletConditions()
            print mlHamiltonJacobi.modelList[l].dirichletGlobalNodeSet
            mlHamiltonJacobi.uList[l][:] = mlHamiltonJacobi.modelList[l].u.dof
            mlHamiltonJacobi.modelList[l].updateCoefficients()
            mlHamiltonJacobi.modelList[l].coefficients.updateSignFunction(mlHamiltonJacobi.modelList[l].q['u'].flat,
                                                                          Numeric.reshape(mlHamiltonJacobi.modelList[l].q['grad(u)'].flat,
                                                                                          (mlHamiltonJacobi.modelList[l].nQuadraturePoints_global,
                                                                                           mlHamiltonJacobi.modelList[l].nSpace_global)),
                                                                          1.5*mlHamiltonJacobi.modelList[l].mesh.h)
            mlHamiltonJacobi.modelList[l].updateCoefficients()
        while (tn < T[test]):
            if timeIntegration[test] != NoIntegration:
                mlHamiltonJacobi.chooseDT(DTSET)
                if nSteps == 0:
                    mlHamiltonJacobi.initializeTimeIntegration()
                    mlHamiltonJacobi.initializeTimeIntegration()
                tn += mlHamiltonJacobi.DT
                if tstring != None:
                    sys.stdout.write(eraseTime)
                else:
                    sys.stdout.write('T = %12.5e, tn = ' % T[test])
                tstring='%12.5e' % (tn,)
                sys.stdout.write(tstring)
                sys.stdout.flush()
                nSteps += 1
                testOut = test + ('%4.4i' % nSteps)
            else:
                for l in range(nLevels):
                    mlHamiltonJacobi.modelList[l].updateZeroLevelSetDirichletConditions()
                    print mlHamiltonJacobi.modelList[l].dirichletGlobalNodeSet
                    mlHamiltonJacobi.uList[l][:] = mlHamiltonJacobi.modelList[l].u.dof
                mlHamiltonJacobi.DT = 1.0
                tn=1.0
                nSteps +=1
                testOut = test
            if nonlinearSolverType == 'NLNI':
                nonlinearSolver.solveMultilevel(uList   = 
                                                mlHamiltonJacobi.uList,
                                                rList   = 
                                                mlHamiltonJacobi.rList)
            elif nonlinearSolverType == 'Newton':
                for l in range(nLevels):
                    mlHamiltonJacobi.modelList[l].getResidual(u =
                                                    mlHamiltonJacobi.uList[l],
                                                               r =
                                                    mlHamiltonJacobi.rList[l])
                    nonlinearSolver[l].solve(u = mlHamiltonJacobi.uList[l],
                                             r = mlHamiltonJacobi.rList[l])
            mlHamiltonJacobi.updateTimeHistory()
            #rnorm = wl2Norm(ni.resid ual(),h)
            #print "rnorm"+`rnorm`
            for l in range(nLevels):
                if DG != True:
    #                     nSolPlot.hardcopy(testOut+'_sol.eps', eps=1,enhanced=1,color=1)
    #                      nResPlot = Gnuplot.Gnuplot()
    #                      nResPlot("set terminal x11")
    #                      nResPlot.title(testOut)
    #                      nResPlot.plot(
    #                          Gnuplot.Data(mlHamiltonJacobi.rList[-1],
    #                                       with='linespoints',
    #                                       title='numerical residual'))
        #  #                  nResPlot.hardcopy(testOut+'_res.eps', eps=1,enhanced=1,color=1)
                    if nd[test]==1:
                        solPlot.title(testOut)
                        solPlot.plot(Gnuplot.Data(mlMesh.meshList[l].nodeArray[:,0],
                                                  mlHamiltonJacobi.modelList[l].u.dof,
                                                  with='linespoints',
                                                  title='numerical solution'))
                        if analyticalSolution[test] != None:
                            aSolPlot.title(testOut)
                            aSolPlot.plot(Gnuplot.Data(Px,
                                                       [analyticalSolution[test].uOfXT(x,tn) for x in P],
                                                       with='lines',
                                                       title='analytical solution'))
                    elif nd[test]==2:
                        nSol = Numeric.reshape(mlHamiltonJacobi.modelList[l].u.dof,(nx[l],ny[l]))
                        solPlot('set parametric')
                        solPlot('set data style lines')
                        solPlot('set hidden')
                        solPlot('set contour base')
                        #solPlot('set cntrparam levels incremental 0.0,0.1,1.0')
                        solPlot.xlabel('x')
                        solPlot.ylabel('y')
                        solPlot.splot(Gnuplot.GridData(nSol,
                                                       x[l],
                                                       y[l],
                                                       binary=0,
                                                       inline=0,
                                                       ))
                        if analyticalSolution[test] != None:
                            for i in range(nx[l]):
                                for j in range(ny[l]):
                                    aSol[l][i,j] = analyticalSolution[test].uOfXT(Numeric.array([x[l][i],y[l][j],0.0]),tn)
                            aSolPlot('set parametric')
                            aSolPlot('set data style lines')
                            aSolPlot('set hidden')
                            aSolPlot('set contour base')
                            #aSolPlot('set cntrparam levels incremental 0.0,0.1,1.0')
                            aSolPlot.xlabel('x')
                            aSolPlot.ylabel('y')
                            aSolPlot.splot(Gnuplot.GridData(aSol[l],
                                                           x[l],
                                                           y[l],
                                                           binary=0,
                                                           inline=0
                                                            ))
            if computeEigenvalues:
                gevals("set terminal x11")
                gevals.title(testOut+' eigenvalues');
                gevals.xlabel(r'real(\lambda)')
                gevals.ylabel(r'imag(\lambda)')
                gevals.plot()
                kplot("set terminal x11")            
                kplot.plot()
                kplot.title(testOut+r' ln(\kappa) vs. ln(1/h)')
                kplot.xlabel(r'ln(\kappa)')
                kplot.xlabel('ln(1/h)')
                kList=[]
            if computeSpaceTimeError or tn >= T[test]:
                eCoarse=1.0
                eFine=1.0
                hCoarse=1.0
                hFine=1.0
                mFine = mlHamiltonJacobi.modelList[-1]
                for l in range(nLevels):
                    utmp[l].dof[:]=0.0
                if analyticalSolution[test] != None:
                    for eN in range(mFine.q['x'].shape[0]):
                        for k in range(mFine.q['x'].shape[1]):
                            asolq[eN,k] = analyticalSolution[test].uOfXT(mFine.q['x'][eN,k],tn)
                for m,jac,mesh,l in zip(mlHamiltonJacobi.modelList,
                                        mlHamiltonJacobi.jacobianList,
                                        mlMesh.meshList,
                                        range(nLevels)):
                    utmp[l].dof[:] = m.u.dof
                    if l < nLevels-1:
                        for lf in range(l,nLevels-1):
                            mlHamiltonJacobi.meshTransfers.prolong_bcList[lf+1].matvec(utmp[lf].dof,utmp[lf+1].dof)
                            #Load the Dirichlet conditions into the projection 
                            for dofN,g in mlHamiltonJacobi.modelList[lf+1].dirichletConditions.DOFBoundaryConditionDict.iteritems():
                                utmp[lf+1].dof[dofN] = g(mlHamiltonJacobi.modelList[lf+1].dirichletConditions.DOFBoundaryPointDict[dofN],tn)
                        #solPlot.replot(Gnuplot.Data(mlMesh.meshList[-1].nodeArray[:,0],
                        #                          utmp[-1].dof,
                        #                          with='linespoints',
                        #                          title='numerical solution'))
                        utmp[-1].getValues(mFine.q['v'],uqtmp)
                    else:
                        uqtmp[:]=mFine.q['u']
                    if analyticalSolution[test] == None:
                        asolq[:]=mFine.q['u']
                    eCoarse=eFine
                    hCoarse=hFine
                    hFine = mesh.h
                    eFine = L2errorSFEM(mFine.q['dx_h'],asolq,uqtmp)
                    if eSpaceTime.has_key(hFine):
                        eSpaceTime[hFine] += mlHamiltonJacobi.DT*eFine**2
                    else:
                        eSpaceTime[hFine] = mlHamiltonJacobi.DT*eFine**2
                    eSpace[hFine] = eFine
                if computeEigenvalues:
                    dev = DenseEigenvalues(jac)
                    dev.computeEigenvalues()
                    try:
                        gevals.replot(
                            Gnuplot.Data([l for l in dev.eigenvalues],
                                         [0.0 for l in dev.eigenvalues],
                                         title='h = %12.5e' % hFine,with='points'))
                    except TypeError:
                        gevals.replot(
                            Gnuplot.Data([l.real for l in dev.eigenvalues],
                                         [l.imag for l in dev.eigenvalues],
                                         title='h = %12.5e' % hFine,
                                         with='points'))
                    k = max(abs(dev.eigenvalues))/min(abs(dev.eigenvalues))
                    kList.append(k)
                    ratio = k*(hFine**2)
                    print "k*h**2 %12.5E" % ratio
        hFine = 0
        errors='$\|error \|_{L_2(\Omega)}$'
        errorsSpaceTime='$\| \|error \|_{L_2 (\Omega) } \|_{L_2(T)}$ '
        orders='spatial order'
        columns='c|'
        hs="$\\Delta t = %4.2e$, $h=$" % mlHamiltonJacobi.DT
        for mesh in mlMesh.meshList:
            columns+='c'
            hCoarse=hFine
            hFine = mesh.h
            hs += "& %4.2e" % hFine
            if hCoarse != 0:
                if eSpace[hFine] != 0.0 and eSpace[hCoarse] != 0.0:
                    p = (log(eSpace[hFine]) - log(eSpace[hCoarse]))/(log(hFine) - log(hCoarse))
                else:
                    p=0
            else:
                p = 0
            errors+="& %4.2e" % eSpace[hFine]
            orders+="& %4.2e" % p
            errorsSpaceTime+="& %4.2e" % sqrt(eSpaceTime[hFine])
        hs += '\\\ \n \hline'
        errors += '\\\ \n'
        orders += '\\\ \n'
        errorsSpaceTime += '\\\ \n'
        report.write('\\begin{center} \n \\begin{tabular}{'+columns+'}\n')
        report.write(hs)
        report.write(errors)
        report.write(orders)
        report.write(errorsSpaceTime)
        report.write('\\end{tabular} \\end{center} \n')
        report.flush()
        if computeEigenvalues:
            kplot.replot(
                Gnuplot.Data([log(1.0/mesh.h) for mesh in mlMesh.meshList],
                             [log(k) for k in kList],
                             with='linespoints'))
            gevals.hardcopy(testOut+'_eig.eps',eps=1,enhanced=1,color=1)
            kplot.hardcopy(testOut+'_cond.eps',eps=1,enhanced=1,color=1)
        if nonlinearSolverType == 'NLNI':
            print "**************NLNI Statistics********************"
            print nonlinearSolver.info()
        else:
            for l in range(nLevels):
                print "**************Nonlinear Solver Statistics at level %i********************" % l
                print nonlinearSolver[l].info()                    
        print "nsteps %i" % nSteps
        if DG != True:
            nsolFile = testOut+'_sol.eps'
            asolFile = testOut+'_asol.eps'
            #solPlot.hardcopy(nsolFile, eps=1,enhanced=1,color=1)
            #aSolPlot.hardcopy(asolFile, eps=1,enhanced=1,color=1)
            report.write("\\begin{figure}{%s}\n\epsfig{file=%s,scale=0.65}\epsfig{file=%s,scale=0.65}\n \\end{figure}" % (test,nsolFile,asolFile))
            report.flush()
        mlMesh.meshList[-1].writeMeshEnsight(test,test)
        mlHamiltonJacobi.modelList[-1].u.name='u'
        #mlHamiltonJacobi.modelList[-1].writeBoundaryTermsEnsight(test)
        mlHamiltonJacobi.modelList[-1].u.writeFunctionEnsight(test,append=False)
        raw_input('Please press return to continue... \n')
    report.write('\\end{document}')
    report.flush()
    import os
    #os.spawnlp(os.P_WAIT,'latex','latex','scalarTransportReport.tex')
    #os.spawnlp(os.P_WAIT,'latex','latex','scalarTransportReport.tex')
    #os.spawnlp(os.P_WAIT,'xdvi','xdvi','scalarTransportReport.dvi')

## @}

if __name__ == '__main__':
    runTests()
