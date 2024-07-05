#! /usr/bin/env python
import os

## need to insert more comments
# \ingroup scripts
# \file proteusRunSSO.py
#
# @{
#   \brief driver for multi-model simulations
#

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
    parser.add_option("-w", "--wait",
                      help="stop after each time step",
                      action="store_true",
                      dest="wait",
                      default=False)
    parser.add_option('--probDir',
                      default='.',
                      help="""where to find problem descriptions""")

    (opts,args) = parser.parse_args()
    #modify path to be able to load proteus test problems
    #for now always insert
    probDir = str(opts.probDir)
    if probDir not in sys.path:
        sys.path.insert(0,probDir)
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
    if opts.logLevel > 0:
        Profiling.openLog(pNameAll+".log",opts.logLevel)
    if opts.viewer:
        Viewers.viewerOn(pNameAll,opts.viewer)
    running = True
    while running:
        if opts.interactive:
            userInput = True
            sys.stdout.write("Enter python commands or (s)tart/(q)quit\n>>>")
        elif opts.batchFileName != "":
            userInput = True
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
                profiler.runctx('runRoutines(pNameAll,pNameList,pList,nList,opts)',{'runRoutines':runRoutines},{'pNameAll':pNameAll,'pNameList':pNameList,'pList':pList,'nList':nList,'opts':opts},pNameAll+'_prof')
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
                runRoutines(pNameAll,pNameList,pList,nList,opts)
        os.chdir('../../')
    if opts.viewer:
        input('\nPress return to close windows and exit... \n')
        if Viewers.viewerType == 'matlab':
            Viewers.viewerPipe.write("quit \n")
        #matlab

def runProblems(pNameAll,pNameList,pList,nList,opts):
    Profiling.memory()
    memBase = Profiling.memLast
    elementQuadratureList=[]
    elementBoundaryQuadratureList=[]
    pM = pList[opts.masterModel]
    nM = nList[opts.masterModel]
    for p,n in zip(pList,nList):
        elementQuadratureDict={}
        for I in p.coefficients.elementIntegralKeys:
            elementQuadratureDict[I] = n.elementQuadrature
        if n.subgridError != None:
            for I in p.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',)+I[1:]] = n.elementQuadrature
        if n.shockCapturing != None:
            for ci in n.shockCapturing.components:
                elementQuadratureDict[('numDiff',ci,ci)] = n.elementQuadrature
        if n.massLumping:
            for ci in list(p.coefficients.mass.keys()):
                elementQuadratureDict[('m',ci)] = Quadrature.SimplexLobattoQuadrature(p.nd,1)
            for I in p.coefficients.elementIntegralKeys:
                elementQuadratureDict[('stab',)+I[1:]] = Quadrature.SimplexLobattoQuadrature(p.nd,1)
        elementBoundaryQuadratureDict={}
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

    if nM.timeIntegration == TimeIntegration.NoIntegration:
        n.T = 1.0
    tn = 0.0
    nSteps = 0
    plotOffSet=0
    for p,n,m in zip(pList,nList,mList):
        m.attachModels(mList)
        if p.initialConditions != None:
            m.setInitialConditions(p.initialConditions,tn)
            m.modelList[-1].viewSolution(plotOffSet,': Initial Condition')
            m.modelList[-1].saveSolution()
        m.modelList[-1].timeIntegration.runCFL = n.runCFL
        plotOffSet += m.modelList[-1].coefficients.nc
        if m.modelList[-1].coefficients.vectorComponents != None:
            plotOffSet += 1
    if nM.timeIntegration != TimeIntegration.NoIntegration:
        DTSET = nM.DT
    else:
        DTSET = 1.0
        pM.T = 1.0
    tstring=None
    eraseTime='\b\b\b\b\b\b\b\b\b\b\b\b'
    timeValues = [tn]
    failedFlag=False
    mM = mList[opts.masterModel]
    firstStep=True
    while (tn < pM.T and not failedFlag==True):
        if opts.wait:
            input('\nPress return to continue... \n')
        failedStep = [False for m in mList]
        if nM.timeIntegration != TimeIntegration.NoIntegration:
            mM.chooseDT(DTSET)
            #check to see if tn+dt > T then reset?
            if abs(tn+mM.DT) > abs(pM.T-tn*1.0e-8):
                DTSET=pM.T-tn
                mM.chooseDT(DTSET)
            #endif
            for m in mList:
                m.chooseDT(mM.DT)
                if nSteps == 0:
                    m.initializeTimeIntegration()
            tn += mM.DT
            if tstring != None:
                sys.stdout.write(eraseTime)
            else:
                sys.stdout.write('T = %12.5e, tn = ' % p.T)
            tstring='%12.5e' % (tn,)
            sys.stdout.write(tstring)
            sys.stdout.flush()
            nSteps += 1
            testOut = pNameAll + ('%4.4i' % nSteps)
        else:
            for m in mList:
                m.DT = 1.0
            tn=1.0
            nSteps +=1
            testOut = pNameAll
        plotOffSet=0
        it=0
        for p,n,m,ls,nls,pName in zip(pList,nList,mList,lsList,nlsList,pNameList):
            print("Model = ",pName)
            #1 stage integration by default
            nStages = 1
            if 'nStagesTime' in dir(n):
                nStages = n.nStagesTime
                #print 'setting nStages= ',nStages
            #endif
            for i in range(nStages):
                if len(mList) == 3:
                    if it  == 0 and not firstStep:
                        for mmls,mmrd,uls,urd in zip(mList[1].modelList,mList[2].modelList,mList[1].uList,mList[2].uList):
                            uls.flat[:] = urd.flat
                            mmls.u[0].dof.flat[:] = mmrd.u[0].dof.flat
                            mmls.initializeTimeIntegration()
                    if it  == 2:
                        for mmls,mmrd,uls,urd in zip(mList[1].modelList,mList[2].modelList,mList[1].uList,mList[2].uList):
                            urd.flat[:] = uls.flat
                            mmrd.u[0].dof.flat[:] = mmls.u[0].dof.flat
                            mmrd.initializeTimeIntegration()
                    if pList[it].weakDirichletConditions != None:
                        for mm in mList[it].modelList:
                            mm.dirichletNodeSetList={}
                            mm.dirichletGlobalNodeSet={}
                            mm.dirichletValues={}
                        for ci in pList[it].weakDirichletConditions:
                            for mm in m.modelList:
                                pList[it].weakDirichletConditions[ci](mm)
                firstStep=False
                failedStage=nls.solveMultilevel(uList   =
                                               m.uList,
                                               rList   =
                                               m.rList)
                m.updateStage()
                failedStep[it] = failedStep[it] or failedStage
            it += 1
            failedFlag = True in failedStep
            timeValues.append(tn)
            m.modelList[-1].viewSolution(plotOffSet,': t=%12.5e' % tn)
            m.modelList[-1].saveSolution()
            plotOffSet+=m.modelList[-1].coefficients.nc
            if m.modelList[-1].coefficients.vectorComponents != None:
                plotOffSet += 1
            if opts.ensight:
                for ci in range(p.coefficients.nc):
                    m.modelList[-1].u[ci].femSpace.writeFunctionEnsight(m.modelList[-1].u[ci],pName,append=True)
            m.updateTimeHistory()
            #mwf  not necessary now?
            #if n.conservativeFlux == 'pwc':
            #    m.modelList[-1].getConservationFluxPWC()
            #elif n.conservativeFlux == 'pwl':
            #    m.modelList[-1].getConservationFluxPWL()
            Profiling.logEvent(multilevelLinearSolver.info(),level=5)
            Profiling.logEvent(nls.info(),level=5)
    #for ci in range(p.coefficients.nc):
    #    mlTransport.modelList[-1].u[ci].femSpace.writeFunctionGnuplot(mlTransport.modelList[-1].u[ci],pName+"""-gnu-%d""" % ci)
    #    raw_input('\nPress return to close windows and exit... \n')
if __name__ == '__main__':
    proteusRun(runProblems)

## @}
