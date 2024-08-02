#!/usr/bin/env python
from ScalarTransport import *
from ScalarTransportTests import *
from LevelSetTests import *


"""

test RKDG via ScalarTransport interface with quadrature for simple
advection problems

"""

def buildProblems(testFlag=0,
                  verbose=0):
    """
    build data structures necessary for specifying test problems:

    testFlag says which one to run?
      0 --- LinearAD_DiracIC  (1d)
      1 --- rotating Cone (2d)

    """
    testProblems = []
    nd = {}
    T  = {}
    coefficients = {}
    getInitialConditions   = {}
    getDirichletConditions = {}
    analyticalSolution     = {}
    if testFlag == 1:
        test = 'RotatingCone2D'
        testProblems.append(test)
        nd[test]=2

        getDirichletConditions[test]=getHomogeneousDBC2D
        N=3.0
        analyticalSolution[test] = RotatingCone2D(1.0/8.0)
        #mwf correct one, but DG has problems for some reason
        #coefficients[test]=UnitSquareRotation()
        #mwf works better with this, still have no diffusion
        coefficients[test]=UnitSquareRotationWithDiffusion(A0=0.0)
        T[test]=0.5
        getInitialConditions[test] = analyticalSolution[test]
        coefficients[test].mass = 'linear'
        coefficients[test].advection = 'linear'
        #mwf correct  coefficients[test].diffusion = None
        #mwf correct coefficients[test].potential = None
        #mwf worked better
        coefficients[test].diffusion = 'constant'
        #mwf worked better
        coefficients[test].potential = 'linear'
        coefficients[test].reaction = None
    else:
        #1d linear advection-diffusion with dirac initial condition:
        #
        #u_t + (bu - a u_x)_x = 0; u(0) = 1; u(1) = 0
        #
        test='LinearAD_DiracIC'

        testProblems.append(test)
        nd[test]=1
        getDirichletConditions[test]=getDBC_hom
        #a0=1.0e-4
        a0=1.0e-2
        #a0=1.0
        #a0=0.0
        A0=Numeric.array([[a0]])
        b0=1.0
        #b0=0.0
        B0=Numeric.array([b0])
        C0=1.0
        M0=0.0
        coefficients[test] = LinearADR_ConstantCoefficients(M=1.0,A=A0,B=B0,C=0.0)
        analyticalSolution[test] = LinearAD_DiracIC(b=B0,a=a0,tStart=0.25)
        T[test]=0.1 #0.5
        getInitialConditions[test] = analyticalSolution[test]
        coefficients[test].mass = 'linear'
        coefficients[test].advection = 'linear'
        coefficients[test].diffusion = 'constant'
        coefficients[test].potential = 'linear'
        coefficients[test].reaction = None
    #end else on testFlag
    #just put these in one dictionary so I know what all
    #has to be specified
    problems = {}
    problems['testProblems']       =testProblems
    problems['nd']                 =nd
    problems['T']                  =T
    problems['coefficients']       =coefficients
    problems['initialConditions']  =getInitialConditions
    problems['dirichletConditions']=getDirichletConditions
    problems['analyticalSolution'] =analyticalSolution

    return problems

#end buildProblems

def buildSimParams(test,TimeIntegrationClass,verbose=0):
    """

    define the necessary flags and tolerances for performing test problems

    """

    #mwf debug
    print 'building simulation details for ',test

    computeEigenvalues = False
    if computeEigenvalues:
        #set flags for eigenvalue computation
        linearSolverType= levelLinearSolverType = 'DenseLU'
        levelNonlinearSolverType = 'Newton'
        nonlinearSolverType = 'NLNI'
    else:
        linearSolverType = levelLinearSolverType = 'SparseLU'
        levelNonlinearSolverType = 'Newton'
        nonlinearSolverType = 'Newton'
    #end else on evalues

    #tolerances
    tolFac   = 1.0e-4
    linTolFac= 1.0e-2
    #time stepping control
    runCFL   = 0.1
    #order of approximation (1 unless doing SSPRK)
    tOrder   = 2
    #pick finite element spaces
    DG       = True #False
    if DG:
        FemSpace = DG_AffineLinearOnSimplexWithNodalBasis
        conservativeFlux = False
        numericalFlux = True
        stabilization= None
        shockCapturing=None
        #mwf added
        shockCapturingDiffusion = None
        quadratureOrder=3
        preSmooths = None
        postSmooths = None
        cycles = None
        nLevels=1
        if test == 'RotatingCone2D':
            nn  =31
        else:
            nn = 51
    else:
        FemSpace = C0_AffineLinearOnSimplexWithNodalBasis
        conservativeFlux = None#'pwc'
        numericalFlux = None
        stabilization='2'
        shockCapturing= None#'1'
        shockCapturingDiffusion = 0.15
        quadratureOrder=3
        preSmooths = 2
        postSmooths = 2
        cycles = 2
        if test == 'RotatingCone2D':
            nLevels = 3
        else:
            nLevels=6 #1d problem
        nn=3 #number of nodes on the coarsest mesh
    #end if on DG

    #collect run parameters
    par = {}
    par['computeEigenvalues']      =computeEigenvalues
    par['linearSolverType']        =linearSolverType
    par['levelLinearSolverType']   =levelLinearSolverType
    par['levelNonlinearSolverType']=levelNonlinearSolverType
    par['nonlinearSolverType']     =nonlinearSolverType
    par['tolFac']            =tolFac
    par['linTolFac']         =linTolFac
    par['runCFL']            =runCFL
    par['DG']                =DG
    par['FemSpace']          =FemSpace
    par['conservativeFlux']  =conservativeFlux
    par['numericalFlux']     =numericalFlux
    par['stabilization']     =stabilization
    par['shockCapturing']    =shockCapturing
    par['shockCapturingDiffusion']    =shockCapturingDiffusion
    par['quadratureOrder']   =quadratureOrder
    par['preSmooths']      =preSmooths
    par['postSmooths']     =postSmooths
    par['cycles']          =cycles
    par['nLevels']         =nLevels
    par['nn']              =nn
    par['timeIntegration'] =TimeIntegrationClass
    par['fullNewtonFlag']  =False
    par['timeIntOrder']    =tOrder
    #par['']  =
    #

    return par
#end buildSimParams

def buildQuadrature(test,tpars,problems):
    """
    setup numerical quadrature data structures
    """
    quadrature = {}
    gq = SimplexGaussQuadrature(problems['nd'][test])
    gq.setOrder(tpars['quadratureOrder'])
    for integral in OneLevelScalarTransport.integralKeys:
        quadrature[integral] = gq
    #end for
    if tpars['stabilization'] is not None:
        quadrature['stab'] = gq
    if tpars['shockCapturing'] is not None:
        quadrature['numDiff'] = gq
    elementBoundaryQuadrature={}
    ebgq = SimplexGaussQuadrature(problems['nd'][test]-1)
    ebgq.setOrder(tpars['quadratureOrder'])
    for elementBoundaryIntegral in OneLevelScalarTransport.elementBoundaryIntegralKeys:
        elementBoundaryQuadrature[elementBoundaryIntegral] = ebgq
    #end boundary quad integral

    tpars['quadrature']= quadrature
    tpars['elementBoundaryQuadrature']=elementBoundaryQuadrature
    return tpars
#end build quadrature

def buildMultilevelMesh(test,tpars,problems):
    mlMesh = []
    nn = tpars['nn']
    nLevels = tpars['nLevels']
    if problems['nd'][test]==1:
        mlMesh = MultiLevelEdgeMesh(nn,1,1,refinementLevels=nLevels)
    elif  problems['nd'][test]==2:
        mlMesh = MultiLevelTriangularMesh(nn,nn,1,
                                          refinementLevels=nLevels)
    elif problems['nd'][test]==3:
        mlMesh = MultiLevelTetrahedralMesh(nn,nn,nn,
                                           refinementLevels=nLevels)
    #end if on dim

    return mlMesh
#end buildMultilevelMesh

def buildMultiLevelScalarTransport(test,tpars,problems,mlMesh):
    """

    """
    tolList=[]
    linTolList=[]
    for l in range(tpars['nLevels']):
        mlMesh.meshList[l].computeGeometricInfo()
        tolList.append(tpars['tolFac']*(mlMesh.meshList[l].h**2))
        linTolList.append(tpars['linTolFac']*(mlMesh.meshList[l].h**2))
    #end l
    atol = min(tolList)
    lin_atol = min(linTolList)
    if (tpars['computeEigenvalues'] or
        tpars['linearSolverType'] == 'DenseLU'):
        MatType = Mat
        matType = 'dense'
    else:
        MatType = SparseMat
        matType = 'csr'
    #end if
    mlScalarTransport = MultiLevelScalarTransport(
        problems['nd'][test],
        mlMesh,
        tpars['FemSpace'],
        tpars['FemSpace'],
        matType,
        problems['dirichletConditions'][test],
        problems['coefficients'][test],
        tpars['quadrature'],
        tpars['elementBoundaryQuadrature'],
        tpars['stabilization'],
        tpars['shockCapturing'],
        tpars['shockCapturingDiffusion'],
        tpars['conservativeFlux'],
        tpars['numericalFlux'],
        tpars['timeIntegration'],
        tpars['timeIntOrder'])

    tpars['MatType'] =MatType
    tpars['atol']    = atol
    tpars['lin_atol']= lin_atol
    tpars['tolList'] = tolList
    tpars['linTolList']= linTolList

    return mlScalarTransport,tpars
#end build mlScalarTransport


def buildSolvers(test,tpars,problems,mlScalarTransport,verbose=0):
    """
    create linear and nonlinear solvers
    """
    #how loud should nonlinear solver be
    printNLinfo=False
    if verbose > 3:
        printNLinfo=True

    levelLinearSolver = None
    #force linearSolver to be SparseLU
    if tpars['linearSolverType'] != 'SparseLU':
        print 'WARNING setting linearSolverType to SparseLU'
        print 'you need to check MatType to make sure SparseMat'
        tpars['linearSolverType'] = 'SparseLU'
    #end if
    levelLinearSolverList=[]
    for l in range(tpars['nLevels']):
        levelLinearSolverList.append(
            SparseLU(mlScalarTransport.jacobianList[l]))
    #end l
    levelLinearSolver = levelLinearSolverList

    linearSolver = None
    #do just plain Newton
    linearSolver = levelLinearSolver
    for l in range(tpars['nLevels']):
        linearSolver[l].printInfo=False
    #end l
    directSolverFlag=True

    #print "Setting up NonlinearSolver"
    #for levelnonlinear solver to be Newton
    if tpars['levelNonlinearSolverType'] != 'Newton':
        print 'WARNING setting levelNonlinearSolverType to Newton'
        tpars['levelNonlinearSolverType'] = 'Newton'
    #end if
    levelNonlinearSolverList=[]
    for l in range(tpars['nLevels']):
        levelNonlinearSolverList.append(
            Newton(linearSolver=linearSolver[l],
                   F=mlScalarTransport.modelList[l],
                   J=mlScalarTransport.jacobianList[l],
                   rtol_r=tpars['tolList'][l],
                   atol_r=tpars['atol'],
                   maxIts=500,
                   convergenceTest = 'r',
                   printInfo=printNLinfo,
                   fullNewton=tpars['fullNewtonFlag'],
                   directSolver=directSolverFlag))
    #end for l
    #for nonlinear solver to be Newton
    if tpars['nonlinearSolverType'] != 'Newton':
        print 'WARNING setting nonlinearSolverType to Newton!'
        tpars['nonlinearSolverType'] = 'Newton'
    #end if
    nonlinearSolver = levelNonlinearSolverList

    return linearSolver,nonlinearSolver,levelLinearSolver
#end buildSolvers

def computeErrors(eSpace,eSpaceTime,eSpaceLast,
                  tn,mlScalarTransport,mlMesh,
                  test,pars,problems,verbose):
    """
    go through and calculate errors on mesh hierarchy

    """
    eCoarse=1.0
    eFine=1.0
    hCoarse=1.0
    hFine=1.0
    analyticalSolution = problems['analyticalSolution']
    for m,jac,mesh in zip(mlScalarTransport.modelList,
                                  mlScalarTransport.jacobianList,
                                  mlMesh.meshList):
        if analyticalSolution[test] is not None:
            eCoarse=eFine
            hCoarse=hFine
            hFine = mesh.h
            eFine = L2errorSFEMvsAF(analyticalSolution[test],
                                    m.q['x'],
                                    m.q['dx_m'],
                                    m.q['u'],tn)
            if eSpace.has_key(hFine):
                eSpaceLast[hFine] = eSpace[hFine]
                if eSpaceTime.has_key(hFine):
                    eSpaceTime[hFine] +=\
                      mlScalarTransport.DT*0.5*(eSpaceLast[hFine]**2 + eFine**2)
                else:
                    eSpaceTime[hFine] =\
                      mlScalarTransport.DT*0.5*(eSpaceLast[hFine]**2 + eFine**2)
                #end else on spaceTime
            #end if on eSpace
            eSpace[hFine] = eFine
        #end analytical solution not none

    #end for

    if analyticalSolution[test] is not None:
        hFine = 0
        errors='||e||_{2}'
        errorsSpaceTime=''
        orders='|| ||e||_2 ||_2'
        for mesh in mlMesh.meshList:
            hCoarse=hFine
            hFine = mesh.h
            if hCoarse != 0:
                if eSpace[hFine] != 0.0 and eSpace[hCoarse] != 0.0:
                    p = (log(eSpace[hFine]) - log(eSpace[hCoarse]))/(log(hFine) - log(hCoarse))
                else:
                    p=0
            else:
                p = 0
            #end if on hCoarse != 0
            errors+="& %4.2e" % eSpace[hFine]
            orders+="& %4.2e" % p
            if eSpaceTime.has_key(hFine): #mwf added if
                errorsSpaceTime+="& %4.2e" % sqrt(eSpaceTime[hFine])
        #end for
        print errors
        print orders
        print errorsSpaceTime
    #end if analytical solution

    return eSpace,eSpaceTime,eSpaceLast

def plotInitial(tn,test,tpars,problems,
                mlMesh,mlScalarTransport):
    """
    plot initial conditions and analytical solutions
    """
    solPlot = None
    aSolPlot= None
    if tpars['DG'] == False:
        solPlot = Gnuplot.Gnuplot()
        solPlot("set terminal x11")
        aSolPlot = Gnuplot.Gnuplot()
        aSolPlot("set terminal x11")
        if problems['nd'][test] == 1:
            if problems['analyticalSolution'][test] is not None:
                solPlot.title(test)
                nap=101
                dxap=Numeric.array([1.0/(nap - 1.0),0.0,0.0])
                P = [(i*dxap) for i in range(nap)]
                Px = [x[0] for x in P]
                solPlot.plot(Gnuplot.Data(mlMesh.meshList[-1].nodeArray[:,0],
                                          mlScalarTransport.modelList[-1].u.dof,
                                          with='linespoints',
                                          title='numerical solution'),
                             Gnuplot.Data(Px,
                                          [problems['analyticalSolution'][test].uOfXT(x,tn) for x in P],
                                          with='lines',
                                          title='analytical solution'))
                aSolPlot.plot(Gnuplot.Data(Px,
                                           [problems['analyticalSolution'][test].uOfXT(x,tn) for x in P],
                                           with='lines'))
            #end if on analytical solution
        elif problems['nd'][test]==2:
            nx = (tpars['nn']-1)*(2**(tpars['nLevels']-1))+1
            ny = nx
            x = Numeric.arange(nx)/float(nx-1)
            y = Numeric.arange(nx)/float(nx-1)
            nSol = Numeric.reshape(mlScalarTransport.modelList[-1].u.dof,
                                   (nx,ny))
            solPlot('set parametric')
            solPlot('set data style lines')
            solPlot('set hidden')
            solPlot('set contour base')
            solPlot('set cntrparam levels incremental 0.1,0.1,1.0')
            solPlot.xlabel('x')
            solPlot.ylabel('y')
            solPlot.splot(Gnuplot.GridData(nSol,
                                           x,
                                           y,
                                           binary=0,
                                           inline=0))
            if problems['analyticalSolution'][test] is not None:
                aSol = Numeric.zeros((nx,ny),Numeric.Float)
                for i in range(nx):
                    for j in range(ny):
                        aSol[i,j]=problems['analyticalSolution'][test].uOfXT(Numeric.array([x[i],y[j],0.0]),tn)

                aSolPlot('set parametric')
                aSolPlot('set data style lines')
                aSolPlot('set hidden')
                aSolPlot('set contour base')
                aSolPlot('set cntrparam levels incremental 0.1,0.1,1.0')
                aSolPlot.xlabel('x')
                aSolPlot.ylabel('y')
                aSolPlot.splot(Gnuplot.GridData(aSol,
                                               x,
                                               y,
                                               binary=0,
                                               inline=0))
            #end if on analytical solution
        #end if on nd ==2
    #end if on not DG
    return solPlot,aSolPlot

def plotTimeStep(solPlot,aSolPlot,tn,test,tpars,problems,
                 mlMesh,mlScalarTransport,testOut):
    """
    plot initial conditions and analytical solutions
    """
    if solPlot is None or aSolPlot is None:
        return solPlot,aSolPlot
    #end nothing to plot with
    if tpars['DG'] == False:
        if problems['nd'][test] == 1:
            if problems['analyticalSolution'][test] is not None:
                solPlot.title(testOut)
                nap=101
                dxap=Numeric.array([1.0/(nap - 1.0),0.0,0.0])
                P = [(i*dxap) for i in range(nap)]
                Px = [x[0] for x in P]
                solPlot.plot(Gnuplot.Data(mlMesh.meshList[-1].nodeArray[:,0],
                                          mlScalarTransport.modelList[-1].u.dof,
                                          with='linespoints',
                                          title='numerical solution'),
                             Gnuplot.Data(Px,
                                          [problems['analyticalSolution'][test].uOfXT(x,tn) for x in P],
                                          with='lines',
                                          title='analytical solution'))
            else:
                solPlot.title(testOut)
                solPlot.plot(Gnuplot.Data(mlMesh.meshList[-1].nodeArray[:,0],
                                          mlScalarTransport.modelList[-1].u.dof,
                                          with='linespoints',
                                          title='numerical solution'))
            #end if on analytical solution
        elif problems['nd'][test]==2:
            nx = (tpars['nn']-1)*(2**(tpars['nLevels']-1))+1
            ny = nx
            x = Numeric.arange(nx)/float(nx-1)
            y = Numeric.arange(nx)/float(nx-1)
            nSol = Numeric.reshape(mlScalarTransport.modelList[-1].u.dof,
                                   (nx,ny))
            solPlot('set parametric')
            solPlot('set data style lines')
            solPlot('set hidden')
            solPlot('set contour base')
            solPlot('set cntrparam levels incremental 0.1,0.1,1.0')
            solPlot.xlabel('x')
            solPlot.ylabel('y')
            solPlot.splot(Gnuplot.GridData(nSol,
                                           x,
                                           y,
                                           binary=0,
                                           inline=0))
            if problems['analyticalSolution'][test] is not None:
                aSol = Numeric.zeros((nx,ny),Numeric.Float)
                for i in range(nx):
                    for j in range(ny):
                        aSol[i,j]=problems['analyticalSolution'][test].uOfXT(Numeric.array([x[i],y[j],0.0]),tn)
                    #end j
                #end i
                aSolPlot('set parametric')
                aSolPlot('set data style lines')
                aSolPlot('set hidden')
                aSolPlot('set contour base')
                aSolPlot('set cntrparam levels incremental 0.1,0.1,1.0')
                aSolPlot.xlabel('x')
                aSolPlot.ylabel('y')
                aSolPlot.splot(Gnuplot.GridData(aSol,
                                               x,
                                               y,
                                               binary=0,
                                               inline=0))
            #end if on analytical solution
        #end if on nd ==2
    #end if on not DG
    return solPlot,aSolPlot

def plotFinal(solPlot,aSolPlot,tn,test,tpars,problems,
              mlMesh,mlScalarTransport,testOut):
    """
    plot out solution and mesh at last step in a couple of formats
    """
    if tpars['DG'] == False:
        solPlot.hardcopy(testOut+'_sol.eps', eps=1,enhanced=1,color=1)
        aSolPlot.hardcopy(testOut+'_asol.eps', eps=1,enhanced=1,color=1)
    #end if
    mlMesh.meshList[-1].writeMeshEnsight(test,test)
    mlScalarTransport.modelList[-1].u.name='u'
    #mlScalarTransport.modelList[-1].writeBoundaryTermsEnsight(test)
    mlScalarTransport.modelList[-1].u.writeFunctionEnsight(test,append=False)

    return solPlot,aSolPlot

if __name__ == '__main__':

    import sys
    import numpy
    from ScalarTransport import *
    from LinearSolvers import *
    from TimeIntegrationTools import *

    verbose = 5
    #testFlag = 0 # LinearAD_Dirac_IC
    testFlag = 1 # rotating clone

    problems = buildProblems(testFlag,verbose)

    test = problems['testProblems'][0] #first test I hope
    #pars = buildSimParams(test,BackwardEuler)
    #pars = buildSimParams(test,ForwardEuler)
    pars = buildSimParams(test,SSPRKintegration)

    pars = buildQuadrature(test,pars,problems)

    mlMesh = buildMultilevelMesh(test,pars,problems)

    mlScalarTransport,pars = buildMultiLevelScalarTransport(test,pars,
                                                            problems,
                                                            mlMesh)

    linearSolver,nonlinearSolver,levelLinearSolver = \
            buildSolvers(test,pars,problems,mlScalarTransport,verbose=verbose)


    #start time loop?

    nstages= pars['timeIntOrder']
    tn = 0.0
    nSteps = 0
    maxSteps= 1000
    eSpace={}
    eSpaceTime={}
    eSpaceLast={}
    mlScalarTransport.setInitialConditions(problems['initialConditions'][test],
                                           tn)
    #end ic set

    solPlot,aSolPlot = plotInitial(tn,test,pars,problems,
                                   mlMesh,mlScalarTransport)


    mlScalarTransport.modelList[-1].timeIntegration.runCFL= pars['runCFL']

    done = False
    while not done:
        mlScalarTransport.chooseDT()

        dtMin = min(problems['T'][test]-tn,mlScalarTransport.DT)
        mlScalarTransport.chooseDT(DTSET=dtMin)
        if nSteps == 0:
            mlScalarTransport.initializeTimeIntegration()
            mlScalarTransport.initializeTimeIntegration()
        #end if
        tn += mlScalarTransport.DT
        print 'taking step to t= ',tn
        nSteps += 1
        testOut = test + ('%4.4i' % nSteps)

        #only Newton iteration for now
        if pars['nonlinearSolverType'] != 'Newton':
            print 'nonlinearSolverType must be Newton'
            sys.exit(1)
        #end if

        #loop through stages
        for s in range(nstages):
            for l in range(pars['nLevels']):
                mlScalarTransport.modelList[l].getResidual(u =
                                                           mlScalarTransport.uList[l],
                                                           r =
                                                           mlScalarTransport.rList[l])
                nonlinearSolver[l].solve(u = mlScalarTransport.uList[l],
                                         r = mlScalarTransport.rList[l])
            #end l loop
            mlScalarTransport.updateStage()

        #end s loop
        print 'max u on fine= ',max(mlScalarTransport.modelList[-1].u.dof.flat)
        print 'min u on fine= ',min(mlScalarTransport.modelList[-1].u.dof.flat)
        mlScalarTransport.modelList[-1].u.name = test
        mlScalarTransport.modelList[-1].u.writeFunctionGnuplot(test,
                                                               append=False)

        if pars['conservativeFlux'] == 'pwc':
            mlScalarTransport.modelList[-1].getConservationFluxPWC()
        elif pars['conservativeFlux'] == 'pwl':
            mlScalarTransport.modelList[-1].getConservationFluxPWL()
        elif pars['numericalFlux'] is not None:
            mlScalarTransport.modelList[-1].e['conservationResidual'].flat[:]=0.0
            for eN in range(mlScalarTransport.modelList[-1].mesh.nElements_global):
                for i in range(mlScalarTransport.modelList[-1].nDOF_element):
                    mlScalarTransport.modelList[-1].e['conservationResidual'][eN]+=mlScalarTransport.modelList[-1].elementResidual[eN,i]
            #end for eN
            #print 'consRes=',mlScalarTransport.modelList[-1].e['conservationResidual']
            print "Max mass cons error "+`max(abs(mlScalarTransport.modelList[-1].e['conservationResidual']))`

        #end numerical flux is not None
        mlScalarTransport.updateTimeHistory()

        solPlot,aSolPlot = plotTimeStep(solPlot,aSolPlot,tn,test,pars,problems,
                                        mlMesh,mlScalarTransport,testOut)

        #compute error
        eSpace,eSpaceTime,eSpaceLast = computeErrors(
            eSpace,eSpaceTime,eSpaceLast,
            tn,mlScalarTransport,mlMesh,test,pars,problems,verbose)

        #figure out if done or not
        done = (abs(tn - problems['T'][test]) < 1.0e-10
                or nSteps >= maxSteps)


    #end while

    solPlot,aSolPlot = plotFinal(solPlot,aSolPlot,tn,test,pars,problems,
                                 mlMesh,mlScalarTransport,testOut)
