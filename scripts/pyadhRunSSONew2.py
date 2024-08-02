#! /usr/bin/env python
import os
import sys
import pickle
import numpy as numpy
from proteus import *
from warnings import *

def proteusRun(runRoutines):
    import optparse
    if sys.version_info[1] >= 5:
        import cProfile as profiler
    else:
        import profile as profiler
    #
    import pstats
    usage = "usage: %prog [options] pFile.py [nFile.py]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-i", "--interactive",
                      help="Read input from stdin",
                      action="store_true",
                      dest="interactive",
                      default='')
    parser.add_option("-M", "--masterModel",
                      help="Set the model that controls the time step",
                      action="store",
                      type="int",
                      dest="masterModel",
                      default=0)
    parser.add_option("-V", "--viewer",
                      help="Set the application to use for viewing results. Can be gnuplot or matlab",
                      action="store",
                      type="string",
                      dest="viewer",
                      default=False)
    parser.add_option("-C", "--plot-coefficients",
                      help="Plot the coefficients of the transport system",
                      action="store_true",
                      dest="plotCoefficients",
                      default=False)
    parser.add_option("-P", "--petsc-options",
                      help="Options for  PETSc",
                      action="store",
                      type="string",
                      dest="petscOptions",
                      default=None)
    parser.add_option("-b", "--batchFile",
                      help="Read input from a file",
                      action="store",
                      type="string",
                      dest="batchFileName",
                      default="")
    parser.add_option("-p", "--profile",
                      help="Generate a profile of the  run",
                      action="store_true",
                      dest="profile",
                      default=False)
    parser.add_option("-m", "--memory",
                      help="Track memory usage of the  run",
                      action="callback",
                      callback=Profiling.memProfOn_callback)
    parser.add_option("-l", "--log",
                      help="Store information about what the code is doing,0=none,10=everything",
                      action="store",
                      type="int",
                      dest="logLevel",
                      default=1)
    parser.add_option("-v", "--verbose",
                      help="Print logging information to standard out",
                      action="callback",
                      callback=Profiling.verboseOn_callback)
    parser.add_option("-E", "--ensight",
                      help="write data in ensight  format",
                      action="store_true",
                      dest="ensight",
                      default=False)
    parser.add_option("-L", "--viewLevels",
                      help="view solution on every level",
                      action="store_true",
                      dest="viewLevels",
                      default=False)
    parser.add_option("-w", "--wait",
                      help="stop after each time step",
                      action="store_true",
                      dest="wait",
                      default=False)
    parser.add_option('--probDir',
                      default='.',
                      help="""where to find problem descriptions""")

    (opts,args) = parser.parse_args()
    if opts.petscOptions != None:
        sys.argv = sys.argv[:-1]+opts.petscOptions.split()
        print(sys.argv)
    comm = Comm.init()
    #mwf modify path to be able to load proteus test problems
    #for now always insert
    probDir = str(opts.probDir)
    if probDir not in sys.path:
        #mwf debug
        #print """inserting probDir= %s into path= %s""" % (probDir,sys.path)
        sys.path.insert(0,probDir)
        #mwf debug
        #print """now path= %s """ % sys.path
    #end if

    if len(args) < 1:
        raise RuntimeError("No input file specified")
    if len(args) > 1:
        raise RuntimeError("Must specify so-file containing p and n filenames")
    pList = []
    nList = []
    pNameList = []
    soFile = open(args[0],'r')
    lines = soFile.readlines()
    for mLine in lines:
        filenames = mLine.split()
        if len(filenames) == 2:
            pList.append(__import__(filenames[0][:-3]))
            pNameList.append(filenames[0][:-3]+filenames[1][:-3])
            nList.append(__import__(filenames[1][:-3]))
    if opts.batchFileName != "":
        inputStream = open(opts.batchFileName,'r')
    else:
        inputStream = sys.stdin
    pNameAll = ''
    for pName in pNameList:
        pNameAll+=pName
    pNameAll=pNameAll[:20]
    if opts.logLevel > 0:
        Profiling.openLog(pNameAll+".log",opts.logLevel)
    if opts.viewer:
        Viewers.viewerOn(pNameAll,opts.viewer)
    #mwf now try to insert a set of simulation flags that may or may not be identical for all models
    simFlagsList = []
    for i in range(len(pNameList)):
        simFlagsList.append({})
    for pName,p,simFlags in zip(pNameList,pList,simFlagsList):
        simFlags['simulationName'] = pName
        simFlags['dataFile']       = pName+'.dat'
        #components to track
        simFlags['components'] = [ci for ci in range(p.coefficients.nc)]
        #visualization section
        if opts.viewer != False:
            simFlags['plotTimes'] = ['All'] #plot solution at each macro time level (also 'Last')
            simFlags['plotQuantities'] = ['u']  #plot soln and exact soln if exists
            p.coefficients.plotCoefficients = opts.plotCoefficients
        #viewers section
    #pName,simFlags defaults
    #figure out how to do batch file,
    #just have to write batch file in terms of simFlagsList[modelNumber] instead of simFlags?
    if opts.batchFileName != "":
        #try to split batch file into executable blocks delimited by
        #s or start on line by itself
        #terminate batch file with a q or quit on a line by itself
        batchfile = inputStream.readlines()
        batchBlocks = []; curBlock = 0
        batchBlocks.append([])
        curline = 0; done = False
        while curline < len(batchfile) and not done:
            line = batchfile[curline]
            newBlock = False
            if not line.isspace():
                newBlock = (line.split()[0] == 's' or
                            line.split()[0] == 'start')
            if newBlock:
                batchBlocks[-1].append('#end of block\n')
                batchBlocks.append([])
            if not line.isspace():
                done = (line.split()[0] == 'q' or
                        line.split()[0] == 'quit')
            if not done and not newBlock:
                batchBlocks[-1].append(line)
            #
            curline += 1
        #end while
        #always pop the last block so have to have a start at the end
        batchBlocks.pop()

        #mwf debug
        #print """end of batch read size batchBlocks= %d """ % len(batchBlocks)
        #for i,block in enumerate(batchBlocks):
        #    print """block[%d]= %s """ % (i,block)

    #end batch
    running = True

    while running:
        if opts.interactive:
            userInput = True
            sys.stdout.write("Enter python commands or (s)tart/(q)quit\n>>>")
        elif opts.batchFileName != "":
            userInput = False
            lines = '\n'.join(batchBlocks.pop(0))
            exec(lines)
            run     = True
            running = len(batchBlocks) > 0
        else:
            userInput = False
            run = True
            running = False
        while (userInput):
            line = inputStream.readline()
            if line:
                if (line.split()[0] == 's' or
                    line.split()[0] == 'start'):
                    userInput = False
                    run = True
                elif (line.split()[0] == 'q' or
                      line.split()[0] == 'quit'):
                    userInput = False
                    run = False
                    running = False
                else:
                    userInput = True
                    exec(line)
                    sys.stdout.write(">>>")
            else:
                userInput = False
                run = False
                running = False
        if run:
            if opts.profile:
                profiler.runctx('runRoutines(pNameAll,pNameList,pList,nList,opts,simFlagsList)',{'runRoutines':runRoutines},{'pNameAll':pNameAll,'pNameList':pNameList,'pList':pList,'nList':nList,'opts':opts,'simFlagsList':simFlagsList},pNameAll+'_prof')
                stats = pstats.Stats(pNameAll+'_prof')
                stats.strip_dirs()
                stats.dump_stats(pNameAll+'_prof_c')
                stats.sort_stats('cumulative')
                if Profiling.verbose:
                    stats.print_stats(30)
                stats.sort_stats('time')
                if Profiling.verbose:
                    stats.print_stats(30)
            else:
                runRoutines(pNameAll,pNameList,pList,nList,opts,simFlagsList)
        os.chdir('../../')
    if opts.viewer:
        input('\nPress return to close windows and exit... \n')
        if Viewers.viewerType == 'matlab':
            Viewers.viewerPipe.write("quit \n")
        #matlab

def runProblems(pNameAll,pNameList,pList,nList,opts,simFlagsList=None):
    Profiling.memory()
    memBase = Profiling.memLast
    elementQuadratureList=[]
    elementBoundaryQuadratureList=[]
    pM = pList[opts.masterModel]
    nM = nList[opts.masterModel]
    if opts.viewer:
        for p in pList:
            p.coefficients.plotCoefficients=opts.plotCoefficients
    for p,n in zip(pList,nList):
        elementQuadratureDict={}
        elemQuadIsDict = isinstance(n.elementQuadrature,dict)
        if elemQuadIsDict: #set terms manually
            for I in p.coefficients.elementIntegralKeys:
                if I in n.elementQuadrature:
                    elementQuadratureDict[I] = n.elementQuadrature[I]
                else:
                    elementQuadratureDict[I] = n.elementQuadrature['default']
        else:
            for I in p.coefficients.elementIntegralKeys:
                elementQuadratureDict[I] = n.elementQuadrature
        if n.subgridError != None:
            for I in p.coefficients.elementIntegralKeys:
                if elemQuadIsDict:
                    if I in n.elementQuadrature:
                        elementQuadratureDict[('stab',)+I[1:]] = n.elementQuadrature[I]
                    else:
                        elementQuadratureDict[('stab',)+I[1:]] = n.elementQuadrature['default']
                else:
                    elementQuadratureDict[('stab',)+I[1:]] = n.elementQuadrature
        if n.shockCapturing != None:
            for ci in n.shockCapturing.components:
                if elemQuadIsDict:
                    if ('numDiff',ci,ci) in n.elementQuadrature:
                        elementQuadratureDict[('numDiff',ci,ci)] = n.elementQuadrature[('numDiff',ci,ci)]
                    else:
                        elementQuadratureDict[('numDiff',ci,ci)] = n.elementQuadrature['default']
                else:
                    elementQuadratureDict[('numDiff',ci,ci)] = n.elementQuadrature
        if n.massLumping:
            for ci in list(p.coefficients.mass.keys()):
                elementQuadratureDict[('m',ci)] = Quadrature.SimplexLobattoQuadrature(p.nd,1)
            for I in p.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',)+I[1:]] = Quadrature.SimplexLobattoQuadrature(p.nd,1)
        if n.reactionLumping:
            for ci in list(p.coefficients.mass.keys()):
                elementQuadratureDict[('r',ci)] = Quadrature.SimplexLobattoQuadrature(p.nd,1)
            for I in p.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',)+I[1:]] = Quadrature.SimplexLobattoQuadrature(p.nd,1)
        elementBoundaryQuadratureDict={}
        if isinstance(n.elementBoundaryQuadrature,dict): #set terms manually
            for I in p.coefficients.elementBoundaryIntegralKeys:
                if I in n.elementBoundaryQuadrature:
                    elementBoundaryQuadratureDict[I] = n.elementBoundaryQuadrature[I]
                else:
                    elementBoundaryQuadratureDict[I] = n.elementBoundaryQuadrature['default']
        else:
            for I in p.coefficients.elementBoundaryIntegralKeys:
                elementBoundaryQuadratureDict[I] = n.elementBoundaryQuadrature
        elementQuadratureList.append(elementQuadratureDict)
        elementBoundaryQuadratureList.append(elementBoundaryQuadratureDict)
    Profiling.logEvent("Setting up MultilevelMesh",level=1)
    mlMesh = None
    mlMeshFileName = pNameAll+"_mesh%dD.%d" % (pM.nd,nM.nLevels)
#     try:
#         mlMeshFile = open(mlMeshFileName,'rb')
#         Profiling.logEvent("Reading mesh",level=2)
#         mlMesh = cPickle.load(mlMeshFile)
#     except:
    if True:
        for n in nList:
            n.nn = nM.nn
            n.nLevels = nM.nLevels
        Profiling.logEvent("Generating mesh",level=2)
        mlMeshFile = open(mlMeshFileName,'wb')
        if pM.nd==1:
            mlMesh = MeshTools.MultilevelEdgeMesh(nM.nn,1,1,
                                                  pM.L[0],pM.L[1],pM.L[2],
                                                  refinementLevels=n.nLevels)
        elif pM.nd==2:
            if pM.polyfile != None:
                print("reading triangle mesh")
                tmesh = TriangleTools.TriangleBaseMesh(baseFlags="q30Den",
                                                       nbase=1,
                                                       verbose=10)
                tmesh.readFromPolyFile(pM.polyfile)
                mesh=tmesh.convertToPyadhMesh(verbose=1)
                mlMesh = MeshTools.MultilevelTriangularMesh(2,2,1)
                mlMesh.meshList[0] = mesh
                #mlMesh.meshList[-1].writeEdgesGnuplot('refinedTriangleMesh')
                #mlMesh.meshList[-1].viewMeshGnuplot('refinedTriangleMesh')
                for l in range(1,nM.nLevels):
                    #print "hello"
                    mlMesh.refine()
#                 print mlMesh.meshList[-1].nodeArray
#                 print mlMesh.meshList[-1].elementNodesArray
#                mlMesh.meshList[-1].writeEdgesGnuplot('refinedTriangleMesh')
#                mlMesh.meshList[-1].viewMeshGnuplot('refinedTriangleMesh')
            else:
                mlMesh = MeshTools.MultilevelTriangularMesh(nM.nn,nM.nn,1,
                                                            pM.L[0],pM.L[1],pM.L[2],
                                                            refinementLevels=nM.nLevels)
        elif pM.nd==3:
            mlMesh = MeshTools.MultilevelTetrahedralMesh(nM.nn,nM.nn,nM.nn,
                                                         pM.L[0],pM.L[1],pM.L[2],
                                                         refinementLevels=nM.nLevels)
        pickle.dump(mlMesh,mlMeshFile,protocol=pickle.HIGHEST_PROTOCOL)
    #mwf debug
    for l in range(n.nLevels):
        print("""proteusRun debug: mesh level=%d  nElements= %d nNodes=%d nFaces=%d diameter= %g """ % (l,
                                                                                      mlMesh.meshList[l].nElements_global,mlMesh.meshList[l].nNodes_global,mlMesh.meshList[l].nElementBoundaries_global,mlMesh.meshList[l].h))


    Profiling.memory("mesh")
    Profiling.logEvent("Setting up MultilevelTransport",level=1)
    tolList=[]
    linTolList=[]
    for l in range(nM.nLevels):
        mlMesh.meshList[l].computeGeometricInfo()
        tolList.append(nM.tolFac*(mlMesh.meshList[l].h**2))
        linTolList.append(nM.linTolFac*(mlMesh.meshList[l].h**2))
    mList=[]
    lsList=[]
    nlsList=[]
    for mi,p,n in zip(list(range(len(pList))),pList,nList):
        mlTransport = Transport.MultilevelTransport(
            p.nd,
            mlMesh,
            n.femSpaces,
            n.femSpaces,
            n.matrix,
            p.dirichletConditions,
            p.coefficients,
            elementQuadratureList[mi],
            elementBoundaryQuadratureList[mi],
            p.fluxBoundaryConditions,
            p.advectiveFluxBoundaryConditions,
            p.diffusiveFluxBoundaryConditions,
            n.subgridError,
            n.shockCapturing,
            n.conservativeFlux,
            n.numericalFluxType,
            n.timeIntegration)
        mlTransport.modelList[-1].timeIntegration.runCFL = n.runCFL
        mList.append(mlTransport)
        (multilevelLinearSolver,directSolverFlag) = LinearSolvers.multilevelLinearSolverChooser(
            linearOperatorList = mlTransport.jacobianList,
            multilevelLinearSolverType = n.multilevelLinearSolver,
            computeSolverRates=True,
            printSolverInfo=False,
            levelLinearSolverType = n.levelLinearSolver,
            computeLevelSolverRates=True,
            printLevelSolverInfo=False,
            smootherType = n.linearSmoother,
            computeSmootherRates=True,
            printSmootherInfo=False,
            prolongList = mlTransport.meshTransfers.prolongList,
            restrictList = mlTransport.meshTransfers.restrictList,
            connectivityListList = [mlTransport.modelList[l].connectionList for l in range(n.nLevels)],
            relativeToleranceList = linTolList,
            absoluteTolerance = min(linTolList),
            solverMaxIts = 1000,
            cycles=3,
            preSmooths=3,
            postSmooths=3)
        lsList.append(multilevelLinearSolver)
        Profiling.logEvent("Setting up NonlinearSolver",level=2)
        multilevelNonlinearSolver = NonlinearSolvers.multilevelNonlinearSolverChooser(
            nonlinearOperatorList = mlTransport.modelList,
            jacobianList = mlTransport.jacobianList,
            multilevelNonlinearSolverType = n.multilevelNonlinearSolver,
            computeSolverRates=True,
            printSolverInfo=False,
            relativeToleranceList = tolList,
            absoluteTolerance = n.atol,
            levelNonlinearSolverType=n.levelNonlinearSolver,
            computeLevelSolverRates=True,
            printLevelSolverInfo=False,
            smootherType = n.nonlinearSmoother,
            computeSmootherRates=True,
            printSmootherInfo=False,
            preSmooths=3,
            postSmooths=3,
            cycles=3,
            maxSolverIts=n.maxNonlinearIts,
            prolong_bcList = mlTransport.meshTransfers.prolong_bcListDict,
            restrict_bcList = mlTransport.meshTransfers.restrict_bcListDict,
            restrict_bcSumList = mlTransport.meshTransfers.restrict_bcSumListDict,
            prolongList = mlTransport.meshTransfers.prolongList,
            restrictList = mlTransport.meshTransfers.restrictList,
            restrictionRowSumList = mlTransport.meshTransfers.restrictSumList,
            connectionListList=[mlTransport.modelList[l].connectionList for l in range(n.nLevels)],
            linearSolverList=multilevelLinearSolver.solverList,
            linearDirectSolverFlag=directSolverFlag,
            solverFullNewtonFlag=n.fullNewtonFlag,
            levelSolverFullNewtonFlag=n.fullNewtonFlag,
            smootherFullNewtonFlag=n.fullNewtonFlag,
            EWtol=False)
        nlsList.append(multilevelNonlinearSolver)
    Profiling.logEvent(Profiling.memory("Solver memory",className='',memSaved=memBase),level=1)
    Profiling.logEvent("Running Solver",level=1)
    try:
        #os.mkdir('tmp/'+pNameAll)
        os.makedirs('tmp/'+pNameAll)
    except:
        pass
    os.chdir('tmp/'+pNameAll)

    #mwf now try to create a simulation output object for each model
    simOutputList = []
    if simFlagsList != None:
        for p,n,simFlags in zip(pList,nList,simFlagsList):
            simOutputList.append(SimTools.SimulationProcessor(flags=simFlags,nLevels=n.nLevels,
                                                              analyticalSolution=p.analyticalSolution))
        #
    else:
        for i in range(len(pList)):
            simOutputList.append(SimTools.SimulationProcessor())
        #
    #if simflags
    if nM.timeIntegration == TimeIntegration.NoIntegration:
        n.T = 1.0
    tn = 0.0
    nSteps = 0
    tnList = []
    for p in pList:
        tnList.append(tn)
    plotOffSet=0

    #
    timeIntegratorList = []
    pNameM = pNameList[opts.masterModel]
    if opts.ensight:
        mList[0].modelList[-1].u[0].femSpace.writeMeshEnsight(pNameM+'master',pNameM+'master')
    for p,n,m,nls,tn,pName in zip(pList,nList,mList,nlsList,tnList,pNameList):
        m.attachModels(mList)
        if p.initialConditions != None:
            m.setInitialConditions(p.initialConditions,tn)
        #mwf figure out how to handle master model specific output in simOutput
        if opts.ensight:
            m.modelList[-1].u[0].femSpace.writeMeshEnsight(pName,pName)
            for ci in range(p.coefficients.nc):
                m.modelList[-1].u[ci].femSpace.writeFunctionEnsight(m.modelList[-1].u[ci],pName,append=False,firstVariable=(ci==0))
                m.modelList[-1].u[ci].femSpace.writeFunctionHeaderEnsight(m.modelList[-1].u[ci],pName,append=False,firstVariable=( (ci==0) and (pName == pNameList[0]) ),case_filename=pNameM+'master')
        #set weak dirichlet conditions, skip if no initial condition,
        #steady state will get set before first call?
        if p.weakDirichletConditions != None and p.initialConditions != None:
            for mm in m.modelList:
                mm.dirichletNodeSetList={}
                mm.dirichletGlobalNodeSet={}
                mm.dirichletValues={}
            for ci in p.weakDirichletConditions:
                for mm in m.modelList:
                    p.weakDirichletConditions[ci](mm)
            #ci
        #end if on weakDirichletConditions
        timeIntegratorList.append(n.timeIntegrator(m,
                                                   nls,
                                                   n.timeIntegration,
                                                   n))

    #zip
    #basic preprocessing for simulation output now
    for p,n,m,tn,simOutput in zip(pList,nList,mList,tnList,simOutputList):
        simOutput.preprocess(p,n,m,tn)

    tstring=None
    eraseTime='\b\b\b\b\b\b\b\b\b\b\b\b'
    timeValues = [tnList[opts.masterModel]]
    failedFlag=False
    mM = mList[opts.masterModel]
    dt = pList[opts.masterModel].T/nList[opts.masterModel].nDTout
    print("dt",dt)
    firstStep=True
    for n,ti,tn in zip(nList,timeIntegratorList,tnList):
        ti.initialize(n.DT,t0=tn,T=pList[opts.masterModel].T)
        print("proteusRunSSONew2 init stage DT",n.DT)
    for m in mList:
        #mwf debug
        #temporarily turn off if not 1d
        if m.modelList[-1].nSpace_global == 1 and opts.plotCoefficients:
        #if opts.plotCoefficients:
            ctemp=m.modelList[-1].coefficients.allocateDummyCoefficients(c=m.modelList[-1].q)
            m.modelList[-1].coefficients.plotCoefficientFunctions(m.modelList[-1].T,ctemp)
    firstStep = True
    while (tnList[opts.masterModel] < pM.T and not failedFlag==True):
        plotOffSet = 0
        if opts.wait:
            input('\nPress return to continue... \n')
        failedStep = [False for m in mList]

        for it in range(len(mList)):
            print("Model Number =",it)
            for l,mm in enumerate(mList[it].modelList):
                preCopy = mm.coefficients.preStep(tnList[it],firstStep=firstStep)
                if preCopy != None and ('copy_uList') in preCopy and preCopy['copy_uList'] == True:
                    #cek have to use dof here because models might have different strong Dirichlet conditions.
                    #mList[it].uList[l].flat[:] = mList[preCopy['uList_model']].uList[l].flat[:]
                    for u_ci_lhs,u_ci_rhs in zip(list(mList[it].modelList[l].u.values()),list(mList[preCopy['uList_model']].modelList[l].u.values())):
                        u_ci_lhs.dof[:] = u_ci_rhs.dof
                    mList[it].modelList[l].setFreeDOF(mList[it].uList[l])
            if pList[it].weakDirichletConditions != None:
                for mm in mList[it].modelList:
                    mm.dirichletNodeSetList={}
                    mm.dirichletGlobalNodeSet={}
                    mm.dirichletValues={}
                for ci in pList[it].weakDirichletConditions:
                    for mm in mList[it].modelList:
                        pList[it].weakDirichletConditions[ci](mm)
                #ci
            #end if on weakDirichlet conditions
            #mwf hack
            if True:#firstStep == True:
                mList[it].chooseDT(nList[it].DT,tnList[it]+dt)
            #mwf debug
            print("""proteusRunSSONew2 tnMaster=%g tn[it]=%g dt=%g model= %s""" % (tnList[opts.masterModel],tnList[it],dt,it))
            print("tnList[it],dt",tnList[it],dt)
            failedStep[it],tnList[it] = timeIntegratorList[it].calculateSolution(tnList[it],tnList[it]+dt)

            #cek debug
            ##\todo interpolating the fine grid onto coarse grid should be optional. Not always the right thing to do.
            for l in range(nM.nLevels-1,0,-1):
                for cj in range(mList[it].meshTransfers.nc):
                    mList[it].meshTransfers.scaled_restrict_bcListDict[cj][l].matvec(mList[it].modelList[l].u[cj].dof,
                                                                                     mList[it].modelList[l-1].u[cj].dof)
                    mList[it].modelList[l].setFreeDOF(mList[it].uList[l])
                mList[it].modelList[l].calculateCoefficients()
                mList[it].modelList[l].updateTimeHistory(tnList[it]+dt)
            #end cek debug
            for l,mm in enumerate(mList[it].modelList):
                postCopy = mm.coefficients.postStep(tnList[it],firstStep=firstStep)
                if postCopy != None and ('copy_uList') in postCopy and postCopy['copy_uList'] == True:
                    #cek have to use dof here because models might have different strong Dirichlet conditions.
                    #mList[postCopy['uList_model']].uList[l].flat[:] = mList[it].uList[l].flat[:]
                    for u_ci_lhs,u_ci_rhs in zip(list(mList[postCopy['uList_model']].modelList[l].u.values()),list(mList[it].modelList[l].u.values())):
                        u_ci_lhs.dof[:] = u_ci_rhs.dof
                    mList[it].modelList[l].setFreeDOF(mList[it].uList[l])
            #postcopy loop
            plotOffSet = simOutputList[it].processTimeLevel(pList[it],nList[it],mList[it],tnList[it],
                                                            plotOffSet)
        #end model time step loops
        #mwf hack what if don't force use of optDT?
        forceOptDT = False#True
        if forceOptDT:
            for n,ti in enumerate(timeIntegratorList):
                if n != opts.masterModel:
                    ti.DTSET = timeIntegratorList[opts.masterModel].optDT
                    dt = timeIntegratorList[opts.masterModel].optDT
        failedFlag = True in failedStep
        if abs(tnList[opts.masterModel]+dt) > abs(pList[opts.masterModel].T-1.0e-8*tnList[opts.masterModel]):
            dt = pList[opts.masterModel].T-tnList[opts.masterModel]
        if failedFlag:
            print("""proteusRunSSO tn=%g dt=%g failedStep= %s""" % (tnList[opts.masterModel],dt,failedStep))
        if abs(pList[opts.masterModel].T-tnList[opts.masterModel]) < 1.0e-8*tnList[opts.masterModel]:
            tnList[opts.masterModel] = pList[opts.masterModel].T
        for m,p,pName in zip(mList,pList,pNameList):
            if opts.ensight:
                for ci in range(p.coefficients.nc):
                    m.modelList[-1].u[ci].femSpace.writeFunctionEnsight(m.modelList[-1].u[ci],pName,append=True)
        timeValues.append(tnList[opts.masterModel])
        firstStep = False
    #end while
    if opts.ensight:
        for m,pName in zip(mList,pNameList):
            m.modelList[-1].u[0].femSpace.endTimeSeriesEnsight(timeValues,pName,pName)
        m.modelList[-1].u[0].femSpace.endTimeSeriesEnsight(timeValues,pNameM+'master',pNameM+'master')
    for p,n,m,tn,simOutput in zip(pList,nList,mList,tnList,simOutputList):
        simOutput.postprocess(p,n,m,tn)

if __name__ == '__main__':
    proteusRun(runProblems)
