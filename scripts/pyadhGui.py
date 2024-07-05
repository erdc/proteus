#!/usr/bin/env python

##
# \addtogroup scripts
# Simulation scripts and other tools for pre and postprocessing files
#
#
#
# \file proteusRun.py
# @{
#   \ingroup scripts
#   \brief driver for single model simulations
#
#
import os
import sys
import socket
import pickle
import numpy as numpy
print(sys.path)
from proteus import *
from warnings import *
from PyQt4 import QtGui, QtCore

def proteusRun(runRoutine):
    import optparse
    import sys
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
    comm = Comm.init()
    print("is initialized",comm.isInitialized())
    print("done with comm.init")
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
    else:
        p = __import__(args[0][:-3])
        if len(args) == 2:
            pName = args[0][:-3]+args[1][:-3]
            n = __import__(args[1][:-3])
        else:
            pName = args[0][:-3]
            #\todo only need to catch file not found case
            try:
                n = __import__(args[0][:-6]+'n')
            except:
                n = p
        pNameProc = pName + repr(comm.rank())
    if opts.batchFileName != "":
        inputStream = open(opts.batchFileName,'r')
    else:
        inputStream = sys.stdin
    if opts.logLevel > 0:
        Profiling.openLog(pName+".log",opts.logLevel)
    if opts.viewer:
        Viewers.viewerOn(pName,opts.viewer)
    #mwf now try to set some simulation control flags ?
    simFlags = {}
    simFlags['simulationName'] = pName
    simFlags['simulationNameProc']=pNameProc
    #components to track
    simFlags['components'] = [ci for ci in range(p.coefficients.nc)]
    #visualization section
    if opts.viewer != False or opts.ensight:
        simFlags['plotTimes'] = ['All'] #plot solution at each macro time level (also 'Last')
        simFlags['plotQuantities'] = ['u']#,"q:('velocity',0)"]  #plot soln and exact soln if exists
        p.coefficients.plotCoefficients = opts.plotCoefficients
    #
    if opts.ensight:
        if 'plotOptions' in simFlags:
            simFlags['plotOptions']['ensight']['on']=True
        else:
            simFlags['plotOptions'] = {'ensight':{'on':True}}
    #get batch file up front?
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
    #mwf debug
    #running = False
    #cek try to stick a GUI here for building batch files
    def doNothin():
        print("hello")
    class MainWindow(QtGui.QMainWindow):
        def __init__(self):
            QtGui.QMainWindow.__init__(self)
            self.resize(250, 150)
            self.setWindowTitle('toolbar')
            self.exit = QtGui.QAction(QtGui.QIcon('icons/exit.png'), 'Exit', self)
            self.exit.setShortcut('Ctrl+Q')
            self.connect(self.exit, QtCore.SIGNAL('triggered()'), QtCore.SLOT('close()'))
            self.toolbar = self.addToolBar('Exit')
            self.toolbar.addAction(self.exit)
    app = QtGui.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    app.exec_()
    while running:
        if opts.interactive:
            userInput = True
            sys.stdout.write("Enter python commands or (s)tart/(q)quit\n>>>")
        elif opts.batchFileName != "":
            userInput = False
            lines = '\n'.join(batchBlocks.pop(0))
            #exec lines
            run     = True
            running = len(batchBlocks) > 0
        else:
            userInput = False
            run = True
            running = False
        while (userInput):
            line = sys.stdin.readline()
            if line:
                if line.isspace():
                    userInput = True
                elif (line.split()[0] == 's' or
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
                    #exec line
                    sys.stdout.write(">>>")
            else:
                userInput = False
                run = False
                running = False
        if run:
            if opts.profile:
                profiler.runctx('runRoutine(pName,p,n,opts)',{'runRoutine':runRoutine},{'pName':pName,'p':p,'n':n,'opts':opts},pName+'_prof')
                stats = pstats.Stats(pName+'_prof')
                stats.strip_dirs()
                stats.dump_stats(pName+'_prof_c')
                stats.sort_stats('cumulative')
                if Profiling.verbose:
                    stats.print_stats(30)
                stats.sort_stats('time')
                if Profiling.verbose:
                    stats.print_stats(30)
            else:
                runRoutine(pName,p,n,opts,simFlags)
        os.chdir('../../')
    if opts.viewer and (comm.rank()==0):
        input('\nPress return to close windows and exit... \n')
        if Viewers.viewerType == 'matlab':
            Viewers.viewerPipe.write("quit \n")
        #matlab

def runProblem(pName,p,n,opts,simFlags=None):
    from proteus import cmeshTools
    comm = Comm.get()
    pNameProc = pName+repr(comm.rank()) #mwf does this agree with pNameProc above?
    Profiling.memory()


    memBase = Profiling.memLast
    elementQuadratureDict={}
    elemQuadIsDict = isinstance(n.elementQuadrature,dict)
    if elemQuadIsDict: #set terms manually
        for I in p.coefficients.elementIntegralKeys:
            if I in n.elementQuadrature:
                elementQuadratureDict[I] = n.elementQuadrature[I]
            else:
                elementQuadratureDict[I] = n.elementQuadrature['default']
        #mwf include nodal quadrature if n file asks for it?
        if 'NodalQuadrature' in list(n.elementQuadrature.keys()):
            elementQuadratureDict['NodalQuadrature'] = n.elementQuadrature['NodalQuadrature']
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
    Profiling.logEvent("Setting up MultilevelMesh",level=1)
    mlMesh = None
    mlMeshFileName = pName+"_mesh%dD.%d" % (p.nd,n.nLevels)
#     try:
#         mlMeshFile = open(mlMeshFileName,'rb')
#         Profiling.logEvent("Reading mesh",level=2)
#         mlMesh = cPickle.load(mlMeshFile)
#     except:
    if True:
        Profiling.logEvent("Generating mesh",level=2)
        mlMeshFile = open(mlMeshFileName,'wb')
        if p.nd==1:
            mlMesh = MeshTools.MultilevelEdgeMesh(n.nn,1,1,
                                                  p.L[0],p.L[1],p.L[2],
                                                  refinementLevels=n.nLevels)
            #mwf hack
            #import pdb
            #pdb.set_trace()
            #import random
            #tags = [random.randint(0,1) for i in range(mlMesh.meshList[n.nLevels-1].nElements_global)]
            #elementTagArray = numpy.zeros((mlMesh.meshList[n.nLevels-1].nElements_global,),'i')
            #elementTagArray[0] = 1; elementTagArray[-1] = 1;
#             mlMesh.locallyRefine(elementTagArray)
#             n.nLevels += 1
        elif p.nd==2:
            if p.polyfile != None:
                #mwf hack
                #if False:
                #    print "reading triangle mesh node and element files directly"
                #    #could put in system call to triangle for generating from poly
                #    nbase = 1
                #    mesh=MeshTools.TriangularMesh()
                #    mesh.generateFromTriangleFiles(p.polyfile,nbase)
                #    mlMesh = MeshTools.MultilevelTriangularMesh(2,2,1)
                #    mlMesh.generateFromExistingCoarseMesh(mesh,n.nLevels)
                #
                print("reading triangle mesh")
                tmesh = TriangleTools.TriangleBaseMesh(baseFlags=n.triangleOptions,
                                                       nbase=1,
                                                       verbose=10)
                tmesh.readFromPolyFile(p.polyfile)
                mesh=tmesh.convertToPyadhMesh(verbose=1)
                mlMesh = MeshTools.MultilevelTriangularMesh(2,2,1)
                mlMesh.generateFromExistingCoarseMesh(mesh,n.nLevels)
                #previous way of reading in from triangle

                #mlMesh.meshList[-1].writeEdgesGnuplot2('refinedTriangleMesh')
                #mlMesh.meshList[-1].viewMeshGnuplotPipe('refinedTriangleMesh')
            else:
                try:
                    mlMesh = MeshTools.MultilevelTriangularMesh(n.nnx,n.nny,1,
                                                                p.L[0],p.L[1],p.L[2],
                                                                refinementLevels=n.nLevels)
                except:
                    mlMesh = MeshTools.MultilevelTriangularMesh(n.nn,n.nn,1,
                                                                p.L[0],p.L[1],p.L[2],
                                                                refinementLevels=n.nLevels)
        elif p.nd==3:
            if p.polyfile != None:
                if not (os.path.exists(p.polyfile+".ele") and
                        os.path.exists(p.polyfile+".node")):
                    tetcmd = "tetgen -%s %s.poly" % (n.triangleOptions,p.polyfile)
                    #mwf debug
                    print("proteusRun trying to run tetgen with %s " % tetcmd)
                    failed = os.system(tetcmd)
                    elefile = "%s.1.ele" % p.polyfile
                    nodefile = "%s.1.node" % p.polyfile
                    #mwf debug
                    #import pdb
                    #pdb.set_trace()
                    assert os.path.exists(elefile)
                    tmp = "%s.ele" % p.polyfile
                    os.rename(elefile,tmp)
                    assert os.path.exists(tmp)
                    assert os.path.exists(nodefile)
                    tmp = "%s.node" % p.polyfile
                    os.rename(nodefile,tmp)
                    assert os.path.exists(tmp)


                print("reading tetgen mesh node and element files directly")
                nbase = 1
                mesh=MeshTools.TetrahedralMesh()
                mesh.generateFromTetgenFiles(p.polyfile,nbase)
                mlMesh = MeshTools.MultilevelTetrahedralMesh(2,2,2)
                mlMesh.generateFromExistingCoarseMesh(mesh,n.nLevels)
            else:
                print("is initialized",comm.isInitialized())
                mlMesh = MeshTools.MultilevelTetrahedralMesh(n.nn,n.nn,n.nn,
                                                             p.L[0],p.L[1],p.L[2],
                                                             refinementLevels=n.nLevels)
            #mlMesh.meshList[-1].writeEdgesGnuplot2('refinedTetgenMesh')
            #mlMesh.meshList[-1].viewMeshGnuplotPipe('refinedTetgenMesh')
            #mlMesh.meshList[-1].writeTetgenFiles("mesh-last",1)

#         cPickle.dump(mlMesh,mlMeshFile,protocol=cPickle.HIGHEST_PROTOCOL)
    adaptMesh = False#True
    tol = 5.0e-2
    error = tol + 1
    maxLevels = n.nLevels+10
    if adaptMesh:
        elementTagArray = numpy.zeros((mlMesh.meshList[-1].nElements_global,),'i')
        #elementTagArray[mlMesh.meshList[-1].nElements_global/2]=1.0
        elementTagArray.flat[:]=1.0
        #mwf hack
        #elementTagArray.flat[mlMesh.meshList[-1].nElements_global-1]=1.0
        #elementTagArray.flat[21]=1.0
    try:
        #os.mkdir('tmp/'+pName)
        os.makedirs('tmp/'+pName)
    except:
        pass
    #mwf save base directory
    baseDirectory = os.getcwd()
    os.chdir('tmp/'+pName)
    firstLoop = True

    while ((error > tol) and
           (n.nLevels <= maxLevels) and
           (adaptMesh==True or firstLoop == True)):
        firstLoop = False
        if (n.nLevels < maxLevels and adaptMesh ==  True):
            #import random
            #for i in range(mlMesh.meshList[n.nLevels-1].nElements_global):
            #    elementTagArray[i] = random.randint(0,1)
            mlMesh.locallyRefine(elementTagArray)
            #mlMesh.meshList[-1].writeEdgesGnuplot2('locallyRefineMesh')
            #mlMesh.meshList[-1].viewMeshGnuplotPipe('locallyRefineMesh')
            n.nLevels += 1
        else:
            adaptMesh=False
        print("error,tol,error>tol",error,tol,error > tol)
        #mwf debug
        #mlMesh.meshList[-1].writeEdgesGnuplot2('mesh')
        #mlMesh.meshList[-1].viewMeshGnuplotPipe('mesh')
        for l in range(n.nLevels):
            print("""proteusRun debug: mesh level=%d  nElements= %d nNodes=%d nFaces=%d diameter= %g """ % (l,
                                                                                          mlMesh.meshList[l].nElements_global,mlMesh.meshList[l].nNodes_global,mlMesh.meshList[l].nElementBoundaries_global,mlMesh.meshList[l].h))


        Profiling.logEvent(Profiling.memory("mesh","proteusRun"),level=1)
        Profiling.logEvent("Setting up MultilevelTransport",level=1)
        tolList=[]
        linTolList=[]
        pOrder=1
        for ci,space in enumerate(n.femSpaces.values()):
            print("Finite Element Space for comnent ci=",ci,space)
            if (space == FemTools.DG_AffineQuadraticOnSimplexWithNodalBasis or
                space == FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis):
                pOrder = 2
        print("Max Polynomial order",pOrder)
        for l in range(n.nLevels):
            if mlMesh.meshList[l].hasGeometricInfo != True:
                mlMesh.meshList[l].computeGeometricInfo()
            tolList.append(n.tolFac*(mlMesh.meshList[l].h**(pOrder+1)))
            linTolList.append(n.linTolFac*(mlMesh.meshList[l].h**(pOrder+1)))
        mlTransport = Transport.MultilevelTransport(
            p.nd,
            mlMesh,
            n.femSpaces,
            n.femSpaces,
            n.matrix,
            p.dirichletConditions,
            p.coefficients,
            elementQuadratureDict,
            elementBoundaryQuadratureDict,
            p.fluxBoundaryConditions,
            p.advectiveFluxBoundaryConditions,
            p.diffusiveFluxBoundaryConditions,
            n.subgridError,
            n.shockCapturing,
            n.conservativeFlux,
            n.numericalFluxType,
            n.timeIntegration)
        print("Total DOF = ",mlTransport.modelList[-1].nFreeVDOF_global)
        (multilevelLinearSolver,directSolverFlag) = LinearSolvers.multilevelLinearSolverChooser(
            linearOperatorList = mlTransport.jacobianList,
            par_linearOperatorList = mlTransport.par_jacobianList,
            multilevelLinearSolverType = n.multilevelLinearSolver,
            computeSolverRates=True,
            printSolverInfo=False, #mwf debug False,True
            levelLinearSolverType = n.levelLinearSolver,
            computeLevelSolverRates=True,
            printLevelSolverInfo=False,#mwf debug False,
            smootherType = n.linearSmoother,
            computeSmootherRates=True,
            printSmootherInfo=False,
            prolongList = mlTransport.meshTransfers.prolongList,
            restrictList = mlTransport.meshTransfers.restrictList,
            connectivityListList = [mlTransport.modelList[l].sparsityInfo for l in range(n.nLevels)],
            relativeToleranceList = linTolList,
            absoluteTolerance = n.atol,
            solverMaxIts = 2000,
            cycles=n.multigridCycles,
            preSmooths=n.preSmooths,
            postSmooths=n.postSmooths,
            computeEigenvalues=n.computeEigenvalues)
        Profiling.logEvent("Setting up NonlinearSolver",level=2)
        multilevelNonlinearSolver = NonlinearSolvers.multilevelNonlinearSolverChooser(
            mlTransport.modelList,
            mlTransport.jacobianList,
            mlTransport.par_jacobianList,
            duList=mlTransport.duList,
            par_duList=mlTransport.par_duList,
            multilevelNonlinearSolverType = n.multilevelNonlinearSolver,
            computeSolverRates=True,
            printSolverInfo=False, #False,True
            relativeToleranceList = tolList,
            absoluteTolerance = n.atol,
            levelNonlinearSolverType=n.levelNonlinearSolver,
            computeLevelSolverRates=True,
            printLevelSolverInfo=False,
            smootherType = n.nonlinearSmoother,
            computeSmootherRates=True,
            printSmootherInfo=False,
            preSmooths=n.preSmooths,
            postSmooths=n.postSmooths,
            cycles=n.multigridCycles,
            maxSolverIts=n.maxNonlinearIts,
            prolong_bcList = mlTransport.meshTransfers.prolong_bcListDict,
            restrict_bcList = mlTransport.meshTransfers.restrict_bcListDict,
            restrict_bcSumList = mlTransport.meshTransfers.restrict_bcSumListDict,
            prolongList = mlTransport.meshTransfers.prolongList,
            restrictList = mlTransport.meshTransfers.restrictList,
            restrictionRowSumList = mlTransport.meshTransfers.restrictSumList,
            connectionListList=[mlTransport.modelList[l].sparsityInfo for l in range(n.nLevels)],
            linearSolverList=multilevelLinearSolver.solverList,
            linearDirectSolverFlag=directSolverFlag,
            solverFullNewtonFlag=n.fullNewtonFlag,
            levelSolverFullNewtonFlag=n.fullNewtonFlag,
            smootherFullNewtonFlag=n.fullNewtonFlag,
            EWtol=False,
            maxLSits=n.maxLineSearches)
    #     multilevelNonlinearSolver = NonlinearSolvers.multilevelNonlinearSolverChooser(
    #         nonlinearOperatorList = mlTransport.modelList,
    #         jacobianList = mlTransport.jacobianList,
    #         multilevelNonlinearSolverType = n.multilevelNonlinearSolver,
    #         computeSolverRates=True,
    #         printSolverInfo=False,
    #         relativeToleranceList = tolList,
    #         absoluteTolerance = n.atol,
    #         levelNonlinearSolverType=n.levelNonlinearSolver,
    #         computeLevelSolverRates=True,
    #         printLevelSolverInfo=False,
    #         smootherType = n.nonlinearSmoother,
    #         computeSmootherRates=True,
    #         printSmootherInfo=False,
    #         preSmooths=3,
    #         postSmooths=3,
    #         cycles=3,
    #         maxSolverIts=n.maxNonlinearIts,
    #         prolong_bcList = mlTransport.meshTransfers.prolong_bcListDict,
    #         restrict_bcList = mlTransport.meshTransfers.restrict_bcListDict,
    #         restrict_bcSumList = mlTransport.meshTransfers.restrict_bcSumListDict,
    #         prolongList = mlTransport.meshTransfers.prolongList,
    #         restrictList = mlTransport.meshTransfers.restrictList,
    #         restrictionRowSumList = mlTransport.meshTransfers.restrictSumList,
    #         connectionListList=[mlTransport.modelList[l].connectionList for l in range(n.nLevels)],
    #         linearSolverList=multilevelLinearSolver.solverList,
    #         linearDirectSolverFlag=directSolverFlag,
    #         solverFullNewtonFlag=n.fullNewtonFlag,
    #         levelSolverFullNewtonFlag=n.fullNewtonFlag,
    #         smootherFullNewtonFlag=n.fullNewtonFlag,
    #         EWtol=False)
        Profiling.logEvent(Profiling.memory("Solver memory",className='',memSaved=memBase),level=1)
        Profiling.logEvent("Running Solver",level=1)
        if n.timeIntegration == TimeIntegration.NoIntegration:
            n.T = 1.0
        tn = 0.0
        nSteps = 0
        if p.initialConditions != None:
            mlTransport.setInitialConditions(p.initialConditions,tn)
            #mwf ask chris if want to use this or not
            #cek maybe left in ability to view all levels, but got rid of initial conditions plot
            if opts.viewLevels:
                for i in range(len(mlTransport.modelList)):
                    mlTransport.modelList[i].viewSolution(titleModifier=': Initial Condition')
                    mlTransport.modelList[i].saveSolution()
                    if opts.wait:
                        input('\nPress return to close windows and exit... \n')
                #mwf undo plots so not in new window?
                Viewers.windowNumber -= len(mlTransport.modelList)
        if opts.ensight:
            comm = Comm.get()
            if comm.isMaster():
                sosOut = open(pName+'.sos','w')
                header = "FORMAT\n"+"type: master_server gold\n"+"SERVERS\n"+"number of servers: "+repr(comm.size())+"\n\n"
                sosOut.write(header)
                for proc in range(comm.size()):
                    machineEntry = ("machine id: "+socket.gethostname()+"\n"+
                                    "executable: ensight8.server\n"+
                                    "data_path: "+os.getcwd()+"\n"+
                                    "casefile: "+pName+repr(proc)+".case\n\n")
                    sosOut.write(machineEntry)
                sosOut.close()
        mlTransport.modelList[-1].timeIntegration.runCFL = n.runCFL
        #set weak dirichlet conditions
        if p.weakDirichletConditions != None:
            for m in mlTransport.modelList:
                m.dirichletNodeSetList={}
                m.dirichletGlobalNodeSet={}
                m.dirichletValues={}
            for ci in p.weakDirichletConditions:
                for m in mlTransport.modelList:
                    p.weakDirichletConditions[ci](m)
        #end if on weakDirichletConditions
        timeIntegrator = n.timeIntegrator(mlTransport,
                                          multilevelNonlinearSolver,
                                          n.timeIntegration,
                                          n)
        if simFlags != None:
            simOutput = SimTools.SimulationProcessor(flags=simFlags,nLevels=n.nLevels,
                                                     analyticalSolution=p.analyticalSolution)
        else:
            simOutput = SimTools.SimulationProcessor()
        errorEstimator = ErrorEstimators.HierarchicalMeshEstimator(mlTransport)

        simOutput.preprocess(p,n,mlTransport,tn)

        timeIntegrator.initialize(n.DT,t0=tn,T=p.T)
        #mwf debug
        #temporarily turn off if not 1d
        if mlTransport.modelList[-1].nSpace_global == 1 and opts.plotCoefficients==True:
            ctemp=mlTransport.modelList[-1].coefficients.allocateDummyCoefficients(c=mlTransport.modelList[-1].q)
            mlTransport.modelList[-1].coefficients.plotCoefficientFunctions(mlTransport.modelList[-1].T,ctemp)
        timeValues = [tn]
        dt = p.T/n.nDTout
        #mwf debug
        print("""proteusRun T=%g nDTout= %d dt=%g """ % (p.T,n.nDTout,dt))
        failedFlag=False
        while (tn < p.T and not failedFlag==True):
            if opts.wait:
                input('\nPress return to continue... \n')

            #mwf debug
            print("""proteusRun tn=%g dt=%g """ % (tn,dt))
            failedFlag,tn = timeIntegrator.calculateSolution(tn,tn+dt)

            simOutput.processTimeLevel(p,n,mlTransport,tn)
            if adaptMesh:
                error,elementTagArray = errorEstimator.calculate()
            print("error,tol,error>tol",error,tol,error > tol)
            if abs(tn) < abs(p.T) and abs(tn+dt) > abs(p.T - tn*1.0e-8):
                dt= p.T-tn
            assert dt > 1.0e-15, "dt= %s too small " % dt
            #mwf debug
            if failedFlag:
                print("""proteusRun tn=%g dt=%g failed= %s""" % (tn,dt,failedFlag))

            if abs(p.T-tn) < 1.0e-15:
                tn = p.T
           #opts
            timeValues.append(tn)
            Profiling.logEvent(multilevelLinearSolver.info(),level=5)
            Profiling.logEvent(multilevelNonlinearSolver.info(),level=5)
        #end while T
        simOutput.postprocess(p,n,mlTransport,tn)
    #end adapt mesh
    print("nLevels",n.nLevels)
    #mwf hack
    #mlMesh.meshList[-1].writeTriangleFiles("mesh-last",1)
if __name__ == '__main__':
    proteusRun(runProblem)
    Profiling.memorySummary()

## @}
