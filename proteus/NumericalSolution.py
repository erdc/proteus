"""
A hierarchy of classes for managing comlete numerical solution implementations
"""

import os
import numpy
from subprocess import check_call

import LinearSolvers
import NonlinearSolvers
import TriangleTools
import MeshTools
import Profiling
import Transport
import SimTools
import Archiver
import Viewers
from Archiver import ArchiveFlags
import Domain

log = Profiling.logEvent

# Global to control whether the kernel starting is active.
embed_ok = True

class NS_base:  # (HasTraits):
    r"""
    The base class for managing the numerical solution of  PDE's.

    The constructor must build all the objects required by a numerical
    method to approximate the solution over a sequence of time intervals.

    calculateSolution(runName) carries out the numerical solution.

    .. graphviz::
       digraph NumericalSolutionHasA {
       node [shape=record, fontname=Helvetica, fontsize=12];
       NS   [label="NumericalSolution" URL="\ref NumericalSolution", style="filled", fillcolor="gray"];
       mList [label="MultilevelTranportModel [n]" URL="\ref proteus::Transport::MultilevelTransport"];
       nsList [label="NonLinearSolver [n] " URL="\ref proteus::NonLinearSolver"];
       lsList [label="LinearSolver [n] " URL="\ref proteus::LinearSolver"];
       pList [label="Problem Specification [n]" URL="\ref proteus::default_p"];
       nList [label="Numerics Specifiation [n]" URL="\ref proteus::default_n"];
       sList [label="Output Specification [n]" URL="\ref proteus::SimTools"];
       so [label="Coupling Specification " URL="\ref proteus::SO_base"];
       ar [label="Archiver" URL="\ref proteus::AR_base"];
       NS -> pList [arrowhead="normal", style="dashed", color="purple"];
       NS -> nList [arrowhead="normal", style="dashed", color="purple"];
       NS -> so [arrowhead="normal", style="dashed", color="purple"];
       NS -> sList [arrowhead="normal", style="dashed", color="purple"];
       NS -> mList [arrowhead="normal", style="dashed", color="purple"];
       NS -> nsList [arrowhead="normal", style="dashed", color="purple"];
       NS -> lsList [arrowhead="normal", style="dashed", color="purple"];
       NS -> ar [arrowhead="normal", style="dashed", color="purple"];
       }
    """

    def __init__(self, so,pList,nList,sList,opts,simFlagsList=None):
        import Comm
        comm=Comm.get()
        self.comm=comm
        message = "Initializing NumericalSolution for "+so.name+"\n System includes: \n"
        for p in pList:
            message += p.name+"\n"
        log(message)
        self.so=so
        self.pList=pList
        self.nList=nList
        self.opts=opts
        self.simFlagsList=simFlagsList
        self.timeValues={}
        Profiling.memory("Memory used before initializing"+so.name)
        memBase = Profiling.memLast #save current memory usage for later
        if not so.useOneMesh:
            so.useOneArchive=False

        log("Setting Archiver(s)")

        if so.useOneArchive:
            self.femSpaceWritten={}
            tmp  = Archiver.XdmfArchive(opts.dataDir,so.name,useTextArchive=opts.useTextArchive,
                                        gatherAtClose=opts.gatherArchive,hotStart=opts.hotStart)
            self.ar = dict([(i,tmp) for i in range(len(self.pList))])
        elif len(self.pList) == 1:
            self.ar = {0:Archiver.XdmfArchive(opts.dataDir,so.name,useTextArchive=opts.useTextArchive,
                                              gatherAtClose=opts.gatherArchive,hotStart=opts.hotStart)} #reuse so.name if possible
        else:
            self.ar = dict([(i,Archiver.XdmfArchive(opts.dataDir,p.name,useTextArchive=opts.useTextArchive,
                                                    gatherAtClose=opts.gatherArchive,hotStart=opts.hotStart)) for i,p in enumerate(self.pList)])
        #by default do not save quadrature point info
        self.archive_q                 = dict([(i,False) for i in range(len(self.pList))]);
        self.archive_ebq_global        = dict([(i,False) for i in range(len(self.pList))]);
        self.archive_ebqe              = dict([(i,False) for i in range(len(self.pList))]);
        self.archive_pod_residuals = dict([(i,False) for i in range(len(self.pList))]);
        if simFlagsList != None:
            assert len(simFlagsList) == len(self.pList), "len(simFlagsList) = %s should be %s " % (len(simFlagsList),len(self.pList))
            for index in range(len(self.pList)):
                if simFlagsList[index].has_key('storeQuantities'):
                    for quant in filter(lambda a: a != None,simFlagsList[index]['storeQuantities']):
                        recType = quant.split(':')
                        if len(recType) > 1 and recType[0] == 'q':
                            self.archive_q[index] = True
                        elif len(recType) > 1 and recType[0] == 'ebq_global':
                            self.archive_ebq_global[index] = True
                        elif len(recType) > 1 and recType[0] == 'ebqe':
                            self.archive_ebqe[index] = True
                        #
                        elif recType[0] == 'pod_residuals':
                            self.archive_pod_residuals[index]=True
                        else:
                            log("Warning Numerical Solution storeQuantity = %s not recognized won't archive" % quant)
                    #
                #
            #
        #
        log("Setting up MultilevelMesh")
        mlMesh_nList = []
        if so.useOneMesh:
            log("Building one multilevel mesh for all models")
            nListForMeshGeneration=[nList[0]]
            pListForMeshGeneration=[pList[0]]
        else:
            log("Building seperate meshes for each model")
            nListForMeshGeneration=nList
            pListForMeshGeneration=pList

        for p,n in zip(pListForMeshGeneration,nListForMeshGeneration):
            if opts.hotStart:
                p.genMesh = False
                log("Hotstarting, using existing mesh "+p.name)
            else:
                log("Generating mesh for "+p.name)
            #support for old-style domain input
            if p.domain == None:
                if p.nd == 1:
                    p.domain = Domain.RectangularDomain(L=p.L[:1],name=p.name)
                elif p.nd == 2:
                    if p.polyfile != None:
                        p.domain = Domain.PlanarStraightLineGraphDomain(fileprefix=p.polyfile,name=p.polyfile)
                    else:
                        p.domain = Domain.RectangularDomain(L=p.L[:2],name=p.name)
                elif p.nd == 3:
                    if p.polyfile != None:
                        p.domain = Domain.PiecewiseLinearComplexDomain(fileprefix=p.polyfile,name=p.polyfile)
                    elif p.meshfile != None:
                        p.domain = Domain.Mesh3DMDomain(p.meshfile)
                    else:
                        p.domain = Domain.RectangularDomain(L=p.L[:3],name=p.name)
                else:
                    raise RuntimeError("No support for domains in more than three dimensions")
            #now generate meshes, could move to Domain and use polymorphism or MeshTools
            if isinstance(p.domain,Domain.RectangularDomain):
                if p.domain.nd == 1:
                    mlMesh = MeshTools.MultilevelEdgeMesh(n.nn,1,1,
                                                          p.domain.L[0],1,1,
                                                          refinementLevels=n.nLevels,
                                                          nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                          parallelPartitioningType=n.parallelPartitioningType)
                elif p.domain.nd == 2:
                    if (n.nnx == n.nny == None):
                        nnx = nny = n.nn
                    else:
                        nnx = n.nnx
                        nny = n.nny
                    log("Building %i x %i rectangular mesh for %s" % (nnx,nny,p.name))
                    mlMesh = MeshTools.MultilevelTriangularMesh(nnx,nny,1,
                                                                p.domain.L[0],p.domain.L[1],1,
                                                                refinementLevels=n.nLevels,
                                                                nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                                parallelPartitioningType=n.parallelPartitioningType)
                elif p.domain.nd == 3:
                    if (n.nnx == n.nny == n.nnz ==None):
                        nnx = nny = nnz = n.nn
                    else:
                        nnx = n.nnx
                        nny = n.nny
                        nnz = n.nnz
                    log("Building %i x %i x %i rectangular mesh for %s" % (nnx,nny,nnz,p.name))

                    if not hasattr(n,'hex'):
                        n.hex = False

                    if not hasattr(n,'NURBS'):
                        n.NURBS = False

                    if (n.NURBS):
                        mlMesh = MeshTools.MultilevelNURBSMesh(nnx,nny,nnz,
                                                               n.px,n.py,n.pz,
                                                                   p.L[0],p.L[1],p.L[2],
                                                                   refinementLevels=n.nLevels,
                                                                   nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                                   parallelPartitioningType=n.parallelPartitioningType)
                    elif (n.hex):
                        if not hasattr(n,'px'):
                            n.px=0
                            n.py=0
                            n.pz=0
                        mlMesh = MeshTools.MultilevelHexahedralMesh(nnx,nny,nnz,
                                                                   n.px,n.py,n.pz,
                                                                   p.L[0],p.L[1],p.L[2],
                                                                   refinementLevels=n.nLevels,
                                                                   nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                                   parallelPartitioningType=n.parallelPartitioningType)
                    else :
                        mlMesh = MeshTools.MultilevelTetrahedralMesh(nnx,nny,nnz,
                                                                   p.L[0],p.L[1],p.L[2],
                                                                   refinementLevels=n.nLevels,
                                                                   nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                                   parallelPartitioningType=n.parallelPartitioningType)

            elif isinstance(p.domain,Domain.PlanarStraightLineGraphDomain):
                log("Calling Triangle to generate 2D mesh for"+p.name)
                tmesh = TriangleTools.TriangleBaseMesh(baseFlags=n.triangleOptions,
                                                       nbase=1,
                                                       verbose=10)
                if comm.isMaster() and p.genMesh:
                    tmesh.readFromPolyFile(p.domain.polyfile)
                    tmesh.writeToFile(p.domain.polyfile)
                    log("Converting to Proteus Mesh")
                    mesh=tmesh.convertToProteusMesh(verbose=1)
                comm.barrier()
                if not comm.isMaster() or not p.genMesh:
                    mesh = MeshTools.TriangularMesh()
                    mesh.generateFromTriangleFiles(filebase=p.domain.polyfile,base=1)
                mlMesh = MeshTools.MultilevelTriangularMesh(0,0,0,skipInit=True,
                                                            nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                            parallelPartitioningType=n.parallelPartitioningType)
                log("Generating %i-level mesh from coarse Triangle mesh" % (n.nLevels,))
                mlMesh.generateFromExistingCoarseMesh(mesh,n.nLevels,
                                                      nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                      parallelPartitioningType=n.parallelPartitioningType)
            elif isinstance(p.domain,Domain.PiecewiseLinearComplexDomain):
                from subprocess import call
                import sys
                if comm.rank() == 0 and (p.genMesh or not (os.path.exists(p.domain.polyfile+".ele") and
                                                           os.path.exists(p.domain.polyfile+".node") and
                                                           os.path.exists(p.domain.polyfile+".face"))):
                    log("Running tetgen to generate 3D mesh for "+p.name,level=1)
                    tetcmd = "tetgen -%s %s.poly" % (n.triangleOptions,p.domain.polyfile)
                    log("Calling tetgen on rank 0 with command %s" % (tetcmd,))

                    check_call(tetcmd, shell=True)

                    log("Done running tetgen")
                    elefile  = "%s.1.ele" % p.domain.polyfile
                    nodefile = "%s.1.node" % p.domain.polyfile
                    facefile = "%s.1.face" % p.domain.polyfile
                    edgefile = "%s.1.edge" % p.domain.polyfile
                    assert os.path.exists(elefile), "no 1.ele"
                    tmp = "%s.ele" % p.domain.polyfile
                    os.rename(elefile,tmp)
                    assert os.path.exists(tmp), "no .ele"
                    assert os.path.exists(nodefile), "no 1.node"
                    tmp = "%s.node" % p.domain.polyfile
                    os.rename(nodefile,tmp)
                    assert os.path.exists(tmp), "no .node"
                    if os.path.exists(facefile):
                        tmp = "%s.face" % p.domain.polyfile
                        os.rename(facefile,tmp)
                        assert os.path.exists(tmp), "no .face"
                    if os.path.exists(edgefile):
                        tmp = "%s.edge" % p.domain.polyfile
                        os.rename(edgefile,tmp)
                        assert os.path.exists(tmp), "no .edge"
                comm.barrier()
                log("Initializing mesh and MultilevelMesh")
                nbase = 1
                mesh=MeshTools.TetrahedralMesh()
                mlMesh = MeshTools.MultilevelTetrahedralMesh(0,0,0,skipInit=True,
                                                             nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                             parallelPartitioningType=n.parallelPartitioningType)
                if opts.generatePartitionedMeshFromFiles:
                    log("Generating partitioned mesh from Tetgen files")
                    mlMesh.generatePartitionedMeshFromTetgenFiles(p.domain.polyfile,nbase,mesh,n.nLevels,
                                                                  nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                                  parallelPartitioningType=n.parallelPartitioningType)
                else:
                    log("Generating coarse global mesh from Tetgen files")
                    mesh.generateFromTetgenFiles(p.domain.polyfile,nbase,parallel = comm.size() > 1)
                    log("Generating partitioned %i-level mesh from coarse global Tetgen mesh" % (n.nLevels,))
                    mlMesh.generateFromExistingCoarseMesh(mesh,n.nLevels,
                                                          nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                          parallelPartitioningType=n.parallelPartitioningType)
            elif isinstance(p.domain,Domain.MeshTetgenDomain):
                mesh=MeshTools.TetrahedralMesh()
                log("Reading coarse mesh from tetgen file")
                mesh.generateFromTetgenFiles(p.domain.meshfile,1)
                mlMesh = MeshTools.MultilevelTetrahedralMesh(0,0,0,skipInit=True,
                                                             nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                             parallelPartitioningType=n.parallelPartitioningType)
                log("Generating %i-level mesh from coarse Tetgen mesh" % (n.nLevels,))
                mlMesh.generateFromExistingCoarseMesh(mesh,n.nLevels,
                                                      nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                      parallelPartitioningType=n.parallelPartitioningType)
            elif isinstance(p.domain,Domain.Mesh3DMDomain):
                mesh=MeshTools.TetrahedralMesh()
                log("Reading coarse mesh from 3DM file")
                mesh.generateFrom3DMFile(p.domain.meshfile)
                mlMesh = MeshTools.MultilevelTetrahedralMesh(0,0,0,skipInit=True,
                                                             nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                             parallelPartitioningType=n.parallelPartitioningType)
                log("Generating %i-level mesh from coarse 3DM mesh" % (n.nLevels,))
                mlMesh.generateFromExistingCoarseMesh(mesh,n.nLevels,
                                                      nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                      parallelPartitioningType=n.parallelPartitioningType)
            elif isinstance(p.domain,Domain.MeshHexDomain):
                mesh=MeshTools.HexahedralMesh()
                log("Reading coarse mesh from file")
                mesh.generateFromHexFile(p.domain.meshfile)
                mlMesh = MeshTools.MultilevelHexahedralMesh(0,0,0,skipInit=True,
                                                             nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                             parallelPartitioningType=n.parallelPartitioningType)
                log("Generating %i-level mesh from coarse mesh" % (n.nLevels,))
                mlMesh.generateFromExistingCoarseMesh(mesh,n.nLevels,
                                                      nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                      parallelPartitioningType=n.parallelPartitioningType)
            mlMesh_nList.append(mlMesh)
            if opts.viewMesh:
                log("Attempting to visualize mesh")
                try:
                    from proteusGraphical import vtkViewers
                    vtkViewers.ViewMesh(mlMesh.meshList[0],viewMaterialTypes=True)
                    vtkViewers.ViewBoundaryMesh(mlMesh.meshList[0],viewBoundaryMaterialTypes=True)
                except:
                    log("NumericalSolution ViewMesh failed for coarse mesh")
            for l in range(n.nLevels):
                try:
                    log(mlMesh.meshList[l].meshInfo())
                except:
                    log("meshInfo() method not implemented for this mesh type")
                if opts.viewMesh and opts.viewLevels and l > 0:
                    log("Attempting to visualize mesh")
                    try:
                        from proteusGraphical import vtkViewers
                        vtkViewers.ViewMesh(mlMesh.meshList[l],title="mesh level %s " % l,
                                            viewMaterialTypes=True)
                        vtkViewers.ViewBoundaryMesh(mlMesh.meshList[l],title="boundary mesh level %s " % l,
                                                    viewBoundaryMaterialTypes=True)
                    except:
                        log("NumericalSolution ViewMesh failed for mesh level %s" % l)
        if so.useOneMesh:
            for p in pList[1:]: mlMesh_nList.append(mlMesh)
        Profiling.memory("Mesh")
        self.modelList=[]
        self.lsList=[]
        self.nlsList=[]
        self.modelSpinUpList = []
        #
        for p in pList:
            p.coefficients.opts = self.opts
            if p.coefficients.sdInfo == {}:
                for ci,ckDict in p.coefficients.diffusion.iteritems():
                    for ck in ckDict.keys():
                        if not p.coefficients.sdInfo.has_key((ci,ck)):
                            p.coefficients.sdInfo[(ci,ck)] = (numpy.arange(start=0,stop=p.nd**2+1,step=p.nd,dtype='i'),
                                                              numpy.array([range(p.nd) for row in range(p.nd)],dtype='i').flatten())
                            log("Numerical Solution Sparse diffusion information key "+`(ci,ck)`+' = '+`p.coefficients.sdInfo[(ci,ck)]`)

        for p,n,s,mlMesh,index in zip(pList,nList,sList,mlMesh_nList,range(len(pList))):
            if so.needEBQ_GLOBAL:
                n.needEBQ_GLOBAL = True
            if so.needEBQ:
                n.needEBQ = True
            ## \todo clean up tolerances: use rtol_u,atol_u and rtol_res, atol_res; allow scaling by mesh diameter
            ## \todo pass in options = (p,n) instead of using monster ctor signature
            tolList=[]
            linTolList=[]
            for l in range(n.nLevels):
                #if mlMesh.meshList[l].hasGeometricInfo != True:
                #    mlMesh.meshList[l].computeGeometricInfo()

                #fac = (mlMesh.meshList[l].h/mlMesh.meshList[0].h)**2
                fac = 1.0
                tolList.append(n.tolFac*fac)
                linTolList.append(n.linTolFac*fac)

            log("Setting up MultilevelTransport for "+p.name)

            model = Transport.MultilevelTransport(p,n,mlMesh,OneLevelTransportType=p.LevelModelType)
            self.modelList.append(model)
            model.name = p.name
            log("Setting "+model.name+" stepController to "+str(n.stepController))
            model.stepController = n.stepController(model,n)
            Profiling.memory("MultilevelTransport for"+p.name)
            log("Setting up MultilevelLinearSolver for"+p.name)
            #allow options database to set model specific parameters?
            linear_solver_options_prefix = None
            if 'linear_solver_options_prefix' in dir(n):
                linear_solver_options_prefix = n.linear_solver_options_prefix

            (multilevelLinearSolver,directSolverFlag) = LinearSolvers.multilevelLinearSolverChooser(
                linearOperatorList = model.jacobianList,
                par_linearOperatorList = model.par_jacobianList,
                multilevelLinearSolverType = n.multilevelLinearSolver,
                computeSolverRates=n.computeLinearSolverRates,
                printSolverInfo=n.printLinearSolverInfo,
                levelLinearSolverType = n.levelLinearSolver,
                computeLevelSolverRates=n.computeLevelLinearSolverRates,
                printLevelSolverInfo=n.printLevelLinearSolverInfo,
                smootherType = n.linearSmoother,
                computeSmootherRates=n.computeLinearSmootherRates,
                printSmootherInfo=n.printLinearSmootherInfo,
                prolongList = model.meshTransfers.prolongList,
                restrictList = model.meshTransfers.restrictList,
                connectivityListList = [model.levelModelList[l].sparsityInfo for l in range(n.nLevels)],
                relativeToleranceList = linTolList,
                absoluteTolerance = n.l_atol_res,
                solverMaxIts = n.linearSolverMaxIts,
                solverConvergenceTest=n.linearSolverConvergenceTest,
                cycles=n.linearWCycles,
                preSmooths=n.linearPreSmooths,
                postSmooths=n.linearPostSmooths,
                ##\todo logic needs to handle element boundary partition too
                parallelUsesFullOverlap=(n.nLayersOfOverlapForParallel > 0 or n.parallelPartitioningType == MeshTools.MeshParallelPartitioningTypes.node),
                par_duList=model.par_duList,
                solver_options_prefix=linear_solver_options_prefix)
            self.lsList.append(multilevelLinearSolver)
            Profiling.memory("MultilevelLinearSolver for "+p.name)
            log("Setting up MultilevelNonLinearSolver for "+p.name)
            self.nlsList.append(NonlinearSolvers.multilevelNonlinearSolverChooser(
                model.levelModelList,
                model.jacobianList,
                model.par_jacobianList,
                duList=model.duList,
                par_duList=model.par_duList,
                multilevelNonlinearSolverType = n.multilevelNonlinearSolver,
                computeSolverRates=n.computeNonlinearSolverRates,
                solverConvergenceTest=n.nonlinearSolverConvergenceTest,
                levelSolverConvergenceTest=n.levelNonlinearSolverConvergenceTest,
                printSolverInfo=n.printNonlinearSolverInfo,
                relativeToleranceList = tolList,
                absoluteTolerance = n.nl_atol_res,
                levelNonlinearSolverType=n.levelNonlinearSolver,
                computeLevelSolverRates=n.computeNonlinearLevelSolverRates,
                printLevelSolverInfo=n.printNonlinearLevelSolverInfo,
                smootherType = n.nonlinearSmoother,
                computeSmootherRates=n.computeNonlinearSmootherRates,
                printSmootherInfo=n.printNonlinearSmootherInfo,
                preSmooths=n.nonlinearPreSmooths,
                postSmooths=n.nonlinearPostSmooths,
                cycles=n.nonlinearWCycles,
                maxSolverIts=n.maxNonlinearIts,
                prolong_bcList = model.meshTransfers.prolong_bcListDict,
                restrict_bcList = model.meshTransfers.restrict_bcListDict,
                restrict_bcSumList = model.meshTransfers.restrict_bcSumListDict,
                prolongList = model.meshTransfers.prolongList,
                restrictList = model.meshTransfers.restrictList,
                restrictionRowSumList = model.meshTransfers.restrictSumList,
                connectionListList=[model.levelModelList[l].sparsityInfo for l in range(n.nLevels)],
                linearSolverList=multilevelLinearSolver.solverList,
                linearDirectSolverFlag=directSolverFlag,
                solverFullNewtonFlag=n.fullNewtonFlag,
                levelSolverFullNewtonFlag=n.fullNewtonFlag,
                smootherFullNewtonFlag=n.fullNewtonFlag,
                EWtol=n.useEisenstatWalker,
                maxLSits=n.maxLineSearches,
                #\todo need to add logic in multilevel NL solver chooser to account for numerical method's stencil as well
                parallelUsesFullOverlap=(n.nLayersOfOverlapForParallel > 0 or n.parallelPartitioningType == MeshTools.MeshParallelPartitioningTypes.node),
                nonlinearSolverNorm = n.nonlinearSolverNorm))
            model.solver=self.nlsList[-1]
            model.viewer = Viewers.V_base(p,n,s)
            Profiling.memory("MultilevelNonlinearSolver for"+p.name)
            #collect models to be used for spin up
            if index in so.modelSpinUpList:
                self.modelSpinUpList.append(model)
        log("Finished setting up models and solvers")
        if self.opts.save_dof:
            for m in self.modelList:
                for lm in m.levelModelList:
                    for ci in range(lm.coefficients.nc):
                        lm.u[ci].dof_last = lm.u[ci].dof.copy()
        self.archiveFlag= so.archiveFlag
        log("Setting up SimTools for "+p.name)
        self.simOutputList = []
        self.auxiliaryVariables = {}
        if self.simFlagsList != None:
            for p,n,simFlags,model,index in zip(pList,nList,simFlagsList,self.modelList,range(len(pList))):
                self.simOutputList.append(SimTools.SimulationProcessor(flags=simFlags,nLevels=n.nLevels,
                                                                       pFile=p,nFile=n,
                                                                       analyticalSolution=p.analyticalSolution))
                model.simTools = self.simOutputList[-1]
                self.auxiliaryVariables[model.name]= [av.attachModel(model,self.ar[index]) for av in n.auxiliaryVariables]
        else:
            for p,n,s,model,index in zip(pList,nList,sList,self.modelList,range(len(pList))):
                self.simOutputList.append(SimTools.SimulationProcessor(pFile=p,nFile=n))
                model.simTools = self.simOutputList[-1]
                model.viewer = Viewers.V_base(p,n,s)
                self.auxiliaryVariables[model.name]= [av.attachModel(model,self.ar[index]) for av in n.auxiliaryVariables]
        for avList in self.auxiliaryVariables.values():
            for av in avList:
                av.attachAuxiliaryVariables(self.auxiliaryVariables)
        log(Profiling.memory("NumericalSolution memory",className='NumericalSolution',memSaved=memBase))
        if so.tnList == None:
            log("Building tnList from model = "+pList[0].name+" nDTout = "+`nList[0].nDTout`)
            self.tnList=[float(n)*nList[0].T/float(nList[0].nDTout)
                         for n in range(nList[0].nDTout+1)]
        else:
            log("Using tnList from so = "+so.name)
            self.tnList = so.tnList
        log("Time sequence"+`self.tnList`)
        log("Setting "+so.name+" systemStepController to object of type "+str(so.systemStepControllerType))
        self.systemStepController = so.systemStepControllerType(self.modelList,stepExact=so.systemStepExact)
        self.systemStepController.setFromOptions(so)
        log("Finished NumericalSolution initialization")

    ## compute the solution
    def calculateSolution(self,runName):
        log("Setting initial conditions",level=0)
        for index,p,n,m,simOutput in zip(range(len(self.modelList)),self.pList,self.nList,self.modelList,self.simOutputList):
            if self.opts.hotStart:
                log("Setting initial conditions from hot start file for "+p.name)
                tCount = len(self.ar[index].tree.getroot()[-1][-1]) - 1
                time = float(self.ar[index].tree.getroot()[-1][-1][-1][0].attrib['Value'])
                if len(self.ar[index].tree.getroot()[-1][-1]) > 1:
                    dt = time - float(self.ar[index].tree.getroot()[-1][-1][-2][0].attrib['Value'])
                else:
                    log("Only one step in hot start file, setting dt to 1.0")
                    dt = 1.0
                log("Last time step in hot start file was t = "+`time`)
                for lm,lu,lr in zip(m.levelModelList,m.uList,m.rList):
                    for cj in range(lm.coefficients.nc):
                        lm.u[cj].femSpace.readFunctionXdmf(self.ar[index],lm.u[cj],tCount)
                        lm.setFreeDOF(lu)
                        lm.timeIntegration.tLast = time
                        lm.timeIntegration.t = time
                        lm.timeIntegration.dt = dt
                self.tCount = tCount+1
            elif p.initialConditions != None:
                log("Setting initial conditions for "+p.name)
                m.setInitialConditions(p.initialConditions,self.tnList[0])
                #It's only safe to calculate the solution and solution
                #gradients because the models aren't attached yet
                for lm in m.levelModelList:
                    lm.calculateSolutionAtQuadrature()
            else:
                log("No initial conditions provided for model "+p.name)
        if self.opts.hotStart:
            if time >= self.tnList[-1] - 1.0e-5:
                log("Modifying time interval to be tnList[-1] + tnList since tnList hasn't been modified already")
                ndtout = len(self.tnList)
                dtout = (self.tnList[-1] - self.tnList[0])/float(ndtout-1)
                self.tnList = [time + i*dtout for i in range(ndtout)]
                log("New tnList"+`self.tnList`)
            else:
                tnListNew=[time]
                for n,t in enumerate(self.tnList):
                    if time < t-1.0e-8:
                        tnListNew.append(t)
                self.tnList=tnListNew
                log("Hotstarting, new tnList is"+`self.tnList`)
        else:
            self.tCount=0#time step counter
        log("Attaching models and running spin-up step if requested")
        for p,n,m,simOutput in zip(self.pList,self.nList,self.modelList,self.simOutputList):
            m.attachModels(self.modelList)
            if m in self.modelSpinUpList:
                log("Spin-Up Estimating initial time derivative and initializing time history for model "+p.name)
                #now the models are attached so we can calculate the coefficients
                for lm,lu,lr in zip(m.levelModelList,
                                    m.uList,
                                    m.rList):
                    #calculate the coefficients, any explicit-in-time
                    #terms will be wrong
                    lm.getResidual(lu,lr)
                    #post-process velocity
                    #lm.calculateAuxiliaryQuantitiesAfterStep()
                    #load in the initial conditions into time
                    #integration history to get explict terms right
                    lm.initializeTimeHistory()
                    lm.timeIntegration.initializeSpaceHistory()
                    #recalculate coefficients
                    lm.getResidual(lu,lr)
                    #calculate consistent time derivative
                    lm.estimate_mt()
                    #post-process velocity
                    #lm.calculateAuxiliaryQuantitiesAfterStep()
                log("Spin-Up Choosing initial time step for model "+p.name)
                m.stepController.initialize_dt_model(self.tnList[0],self.tnList[1])
                #mwf what if user wants spin-up to be over (t_0,t_1)?

                if m.stepController.stepExact and m.stepController.t_model_last != self.tnList[1]:
                    log("Spin-up step exact called for model %s" % (m.name,),level=3)
                    m.stepController.stepExact_model(self.tnList[1])
                log("Spin-Up Initializing time history for model step controller")
                m.stepController.initializeTimeHistory()
                m.stepController.setInitialGuess(m.uList,m.rList)
                solverFailed = m.solver.solveMultilevel(uList=m.uList,
                                                        rList=m.rList,
                                                        par_uList=m.par_uList,
                                                        par_rList=m.par_rList)
                Profiling.memory("solver.solveMultilevel")
                if solverFailed:
                    log("Spin-Up Step Failed t=%12.5e, dt=%12.5e for model %s, CONTINUING ANYWAY!" %  (m.stepController.t_model,
                                                                                                     m.stepController.dt_model,
                                                                                                     m.name))

                else:
                    if n.restrictFineSolutionToAllMeshes:
                        log("Using interpolant of fine mesh an all meshes")
                        self.restrictFromFineMesh(m)
                    self.systemStepController.modelStepTaken(m,self.tnList[0])
                    log("Spin-Up Step Taken, Model step t=%12.5e, dt=%12.5e for model %s" % (m.stepController.t_model,
                                                                                             m.stepController.dt_model,
                                                                                             m.name))
        for p,n,m,simOutput,index in zip(self.pList,self.nList,self.modelList,self.simOutputList,range(len(self.pList))):
            if not self.opts.hotStart:
                log("Archiving initial conditions")
                self.archiveInitialSolution(m,index)
            else:
                self.ar[index].domain = self.ar[index].tree.find("Domain")
            self.initializeViewSolution(m)
            log("Estimating initial time derivative and initializing time history for model "+p.name)
            #now the models are attached so we can calculate the coefficients
            for lm,lu,lr in zip(m.levelModelList,
                                m.uList,
                                m.rList):
                if self.opts.save_dof:
                    for ci in range(lm.coefficients.nc):
                        lm.u[ci].dof_last[:] = lm.u[ci].dof
                #calculate the coefficients, any explicit terms will be wrong
                lm.timeTerm=False
                lm.getResidual(lu,lr)
                #post-process velocity
                #lm.calculateAuxiliaryQuantitiesAfterStep()
                #load in the initial conditions into time integration history to get explict terms right
                lm.initializeTimeHistory()
                lm.timeIntegration.initializeSpaceHistory()
                #recalculate  coefficients with the explicit terms correct
                lm.getResidual(lu,lr)
                #post-process velocity
                #lm.calculateAuxiliaryQuantitiesAfterStep()
                lm.timeTerm=True
                #calculate consistent
                lm.estimate_mt()
                #
            log("Choosing initial time step for model "+p.name)
            m.stepController.initialize_dt_model(self.tnList[0],self.tnList[1])
            #recalculate  with all terms ready
            for lm,lu,lr in zip(m.levelModelList,
                                m.uList,
                                m.rList):
                lm.getResidual(lu,lr)
            log("Initializing time history for model step controller")
            m.stepController.initializeTimeHistory()
            log("Auxiliary variable calculations for model %s" % (m.name,))
            for av in self.auxiliaryVariables[m.name]:
                av.calculate()
            if not self.opts.cacheArchive:
                self.ar[index].sync()
        self.systemStepController.initialize_dt_system(self.tnList[0],self.tnList[1]) #may reset other dt's
        log("Starting time stepping",level=0)
        self.firstStep = True ##\todo get rid of firstStep flag in NumericalSolution if possible?
        systemStepFailed=False
        stepFailed=False

        #NS_base has a fairly complicated time stepping loop structure
        #to accommodate fairly general split operator approaches. The
        #outer loop is over user defined time intervals for the entire
        #system of models. The next loop is over potentially adaptive
        #steps for the entire system. The next loop is for iterations
        #over the entire system such as for interactive split
        #operator. The next loop is for a sequence of model steps such
        #as for alternating split operator or fractional step
        #schemes. The next loop is for each model to step, potentially
        #adaptively, to the time in the stepSequence. Lastly there is
        #a loop for substeps(stages).

        for (self.tn_last,self.tn) in zip(self.tnList[:-1],self.tnList[1:]):
            log("==============================================================",level=0)
            log("Solving over interval [%12.5e,%12.5e]" % (self.tn_last,self.tn),level=0)
            log("==============================================================",level=0)
            if self.opts.save_dof:
                for m in self.modelList:
                    for lm in m.levelModelList:
                        for ci in range(lm.coefficients.nc):
                            lm.u[ci].dof_last[:] = lm.u[ci].dof
            if self.systemStepController.stepExact and self.systemStepController.t_system_last != self.tn:
                self.systemStepController.stepExact_system(self.tn)
            while self.systemStepController.t_system_last < self.tn:

                log("System time step t=%12.5e, dt=%12.5e" % (self.systemStepController.t_system,
                                                              self.systemStepController.dt_system),level=3)

                while (not self.systemStepController.converged() and
                       not systemStepFailed):
                    log("Split operator iteration %i" % (self.systemStepController.its,),level=3)

                    for (self.t_stepSequence,model) in self.systemStepController.stepSequence:

                        log("Model: %s" % (model.name),level=1)
                        log("Fractional step %12.5e for model %s" % (self.t_stepSequence,model.name),level=3)

                        for m in model.levelModelList:
                            if m.movingDomain and m.tLast_mesh != self.systemStepController.t_system_last:
                                m.t_mesh = self.systemStepController.t_system_last
                                m.updateAfterMeshMotion()
                                m.tLast_mesh = m.t_mesh
                        self.preStep(model)
                        self.setWeakDirichletConditions(model)

                        stepFailed = False
                        if model.stepController.stepExact and model.stepController.t_model_last != self.t_stepSequence:
                            log("Step exact called for model %s" % (model.name,),level=3)
                            model.stepController.stepExact_model(self.t_stepSequence)
                        while (model.stepController.t_model_last < self.t_stepSequence and
                               not stepFailed and
                               not self.systemStepController.exitModelStep[model]):

                            log("Model step t=%12.5e, dt=%12.5e for model %s" % (model.stepController.t_model,
                                                                                 model.stepController.dt_model,
                                                                                 model.name),level=3)

                            for self.tSubstep in model.stepController.substeps:

                                log("Model substep t=%12.5e for model %s" % (self.tSubstep,model.name),level=3)
                                #TODO: model.stepController.substeps doesn't seem to be updated after a solver failure unless model.stepController.stepExact is true
                                log("Model substep t=%12.5e for model %s model.timeIntegration.t= %12.5e" % (self.tSubstep,model.name,model.levelModelList[-1].timeIntegration.t),level=3)

                                model.stepController.setInitialGuess(model.uList,model.rList)

                                solverFailed = model.solver.solveMultilevel(uList=model.uList,
                                                                            rList=model.rList,
                                                                            par_uList=model.par_uList,
                                                                            par_rList=model.par_rList)
                                Profiling.memory("solver.solveMultilevel")
                                if self.opts.wait:
                                    raw_input("Hit any key to continue")
                                if solverFailed:
                                    break
                                else:
                                    if n.restrictFineSolutionToAllMeshes:
                                        log("Using interpolant of fine mesh an all meshes")
                                        self.restrictFromFineMesh(model)
                                    model.stepController.updateSubstep()
                            #end model substeps
                            if solverFailed:
                                log("Step failed due to solver failure")
                                stepFailed = not self.systemStepController.retryModelStep_solverFailure(model)
                            elif model.stepController.errorFailure():
                                log("Step failed due to error failure")
                                stepFailed = not self.systemStepController.retryModelStep_errorFailure(model)
                            else:
                                #set up next step
                                self.systemStepController.modelStepTaken(model,self.t_stepSequence)
                                log("Step Taken, t_stepSequence= %s Model step t=%12.5e, dt=%12.5e for model %s" % (self.t_stepSequence,
                                                                                                                    model.stepController.t_model,
                                                                                                                    model.stepController.dt_model,
                                                                                                                    model.name),level=3)
                                for av in self.auxiliaryVariables[model.name]:
                                    av.calculate()
                        #end model step
                        if stepFailed:
                            log("Sequence step failed")
                            if not self.systemStepController.ignoreSequenceStepFailure(model):
                                break
                            else:
                                log("IGNORING STEP FAILURE")
                                self.postStep(model)
                                self.systemStepController.sequenceStepTaken(model)
                        else:
                            self.postStep(model)
                            self.systemStepController.sequenceStepTaken(model)
                    #end model split operator step
                    if stepFailed:
                        systemStepFailed = not self.systemStepController.retrySequence_modelStepFailure()
                        if not systemStepFailed:
                            stepFailed=False
                            log("Retrying sequence")
                        else:
                            log("Sequence failed")
                    else:
                        self.firstStep=False
                        systemStepFailed=False
                        self.systemStepController.sequenceTaken()
                        for index,model in enumerate(self.modelList):
                            self.viewSolution(model,index)
                        if self.archiveFlag == ArchiveFlags.EVERY_MODEL_STEP:
                            self.tCount+=1
                            for index,model in enumerate(self.modelList):
                                self.archiveSolution(model,index,self.systemStepController.t_system)
                            if not self.opts.cacheArchive:
                                self.ar[index].sync()
                #end system split operator sequence
                if systemStepFailed:
                    log("System Step Failed")
                    #go ahead and update as if the time step had succeeded
                    self.postStep(model)
                    self.systemStepController.modelStepTaken(model,self.t_stepSequence)
                    self.systemStepController.sequenceTaken()
                    self.systemStepController.updateTimeHistory()
                    #you're dead if retrySequence didn't work
                    #go ahead and calculate auxiliary variables for failed solution
                    for model in self.modelList:
                        for av in self.auxiliaryVariables[model.name]:
                            av.calculate()
                    log("Step Failed, Model step t=%12.5e, dt=%12.5e for model %s" % (model.stepController.t_model,
                                                                                      model.stepController.dt_model,
                                                                                      model.name))
                    break
                else:
                    self.systemStepController.updateTimeHistory()
                    self.systemStepController.choose_dt_system()
                    log("Step Taken, System time step t=%12.5e, dt=%12.5e" % (self.systemStepController.t_system,
                                                                              self.systemStepController.dt_system))
                    if self.systemStepController.stepExact and self.systemStepController.t_system_last != self.tn:
                        self.systemStepController.stepExact_system(self.tn)
                    log("Step Taken, Model step t=%12.5e, dt=%12.5e for model %s" % (model.stepController.t_model,
                                                                                     model.stepController.dt_model,
                                                                                     model.name))

                if self.archiveFlag == ArchiveFlags.EVERY_SEQUENCE_STEP:
                    self.tCount+=1
                    for index,model in enumerate(self.modelList):
                        self.archiveSolution(model,index,self.systemStepController.t_system_last)
                    if not self.opts.cacheArchive:
                        self.ar[index].sync()
            #end system step iterations
            if self.archiveFlag == ArchiveFlags.EVERY_USER_STEP:
                self.tCount+=1
                for index,model in enumerate(self.modelList):
                    self.archiveSolution(model,index,self.systemStepController.t_system_last)
                if not self.opts.cacheArchive:
                    self.ar[index].sync()
            if systemStepFailed:
                break
        log("Finished calculating solution",level=3)
        for index,model in enumerate(self.modelList):
            self.finalizeViewSolution(model)
            self.closeArchive(model,index)
        return systemStepFailed
    #
    #try to make preStep and postStep just manipulate "current values" and let the step controllers manage the history setting
    ##intermodel transfer before a solution step
    def preStep(self,model):
        for level,levelModel in enumerate(model.levelModelList):
            preCopy = levelModel.coefficients.preStep(model.stepController.t_model,firstStep=self.firstStep)
            if (preCopy != None and preCopy.has_key(('copy_uList')) and preCopy['copy_uList'] == True):
                for u_ci_lhs,u_ci_rhs in zip(levelModel.u.values(),self.modelList[preCopy['uList_model']].levelModelList[level].u.values()):
                    u_ci_lhs.dof[:] = u_ci_rhs.dof
                levelModel.setFreeDOF(model.uList[level])
            if preCopy != None and preCopy.has_key(('clear_uList')) and preCopy['clear_uList'] == True:
                for u_ci_lhs in levelModel.u.values():
                    u_ci_lhs.dof[:] = 0.0
                levelModel.setFreeDOF(model.uList[level])
            if preCopy != None and preCopy.has_key(('reset_uList')) and preCopy['reset_uList'] == True:
                levelModel.setFreeDOF(model.uList[level])
                levelModel.getResidual(model.uList[level],model.rList[level])

    ##intermodel transfer after a step
    def postStep(self,model):
        for level,levelModel in enumerate(model.levelModelList):
            postCopy = levelModel.coefficients.postStep(model.stepController.t_model,firstStep=self.firstStep)
            if postCopy != None and postCopy.has_key(('copy_uList')) and postCopy['copy_uList'] == True:
                for u_ci_lhs,u_ci_rhs in zip(self.modelList[postCopy['uList_model']].levelModelList[level].u.values(),model.levelModelList[level].u.values()):
                    u_ci_lhs.dof[:] = u_ci_rhs.dof
                self.modelList[postCopy['uList_model']].levelModelList[level].setFreeDOF(self.modelList[postCopy['uList_model']].uList[level])

    def setWeakDirichletConditions(self,model):
        if model.weakDirichletConditions != None:
            for levelModel in model.levelModelList:
                levelModel.dirichletNodeSetList={}
                levelModel.dirichletGlobalNodeSet={}
                levelModel.dirichletValues={}
            for ci in model.weakDirichletConditions:
                for levelModel in model.levelModelList:
                    model.weakDirichletConditions[ci](levelModel)

    def restrictFromFineMesh(self,model):
        for level in range(len(model.levelModelList)-1,0,-1):
            for cj in range(model.levelModelList[-1].coefficients.nc):
                model.meshTransfers.interp_bcListDict[cj][level].matvec(model.levelModelList[level].u[cj].dof,
                                                                                 model.levelModelList[level-1].u[cj].dof)
            model.levelModelList[level-1].setFreeDOF(model.uList[level-1])
            model.levelModelList[level-1].calculateCoefficients()

    ##save model's initial solution values to archive
    def archiveInitialSolution(self,model,index):
        import xml.etree.ElementTree as ElementTree
        if self.archiveFlag == ArchiveFlags.UNDEFINED:
            return
        log("Writing initial mesh for  model = "+model.name,level=3)
        log("Writing initial conditions for  model = "+model.name,level=3)
        if not self.so.useOneArchive or index==0:
            self.ar[index].domain = ElementTree.SubElement(self.ar[index].tree.getroot(),"Domain")
        if self.so.useOneArchive:
            model.levelModelList[-1].archiveFiniteElementSolutions(self.ar[index],self.tnList[0],self.tCount,initialPhase=True,
                                                                   writeVectors=True,meshChanged=True,femSpaceWritten=self.femSpaceWritten,
                                                                   writeVelocityPostProcessor=self.opts.writeVPP)
        else:
            model.levelModelList[-1].archiveFiniteElementSolutions(self.ar[index],self.tnList[0],self.tCount,initialPhase=True,
                                                                   writeVectors=True,meshChanged=True,
                                                                   writeVelocityPostProcessor=self.opts.writeVPP)
        #could just pull the code and flags out from SimTools rathter than asking it to parse them
        #uses values in simFlags['storeQuantities']
        #q dictionary
        if self.archive_q[index] == True:
            scalarKeys = model.simTools.getScalarElementStorageKeys(model,self.tnList[0])
            vectorKeys = model.simTools.getVectorElementStorageKeys(model,self.tnList[0])
            tensorKeys = model.simTools.getTensorElementStorageKeys(model,self.tnList[0])
            model.levelModelList[-1].archiveElementQuadratureValues(self.ar[index],self.tnList[0],self.tCount,
                                                                    scalarKeys=scalarKeys,vectorKeys=vectorKeys,tensorKeys=tensorKeys,
                                                                    initialPhase=True,meshChanged=True)
        if self.archive_ebq_global[index] == True:
            #ebq_global dictionary
            scalarKeys = model.simTools.getScalarElementBoundaryStorageKeys(model,self.tnList[0])
            vectorKeys = model.simTools.getVectorElementBoundaryStorageKeys(model,self.tnList[0])
            tensorKeys = model.simTools.getTensorElementBoundaryStorageKeys(model,self.tnList[0])
            model.levelModelList[-1].archiveElementBoundaryQuadratureValues(self.ar[index],self.tnList[0],self.tCount,
                                                                                scalarKeys=scalarKeys,vectorKeys=vectorKeys,tensorKeys=tensorKeys,
                                                                                initialPhase=True,meshChanged=True)
        if self.archive_ebqe[index] == True:
            #ebqe dictionary
            scalarKeys = model.simTools.getScalarExteriorElementBoundaryStorageKeys(model,self.tnList[0])
            vectorKeys = model.simTools.getVectorExteriorElementBoundaryStorageKeys(model,self.tnList[0])
            tensorKeys = model.simTools.getTensorExteriorElementBoundaryStorageKeys(model,self.tnList[0])
            model.levelModelList[-1].archiveExteriorElementBoundaryQuadratureValues(self.ar[index],self.tnList[0],self.tCount,
                                                                                    scalarKeys=scalarKeys,vectorKeys=vectorKeys,tensorKeys=tensorKeys,
                                                                                    initialPhase=True,meshChanged=True)

        #for nonlinear POD 
        if self.archive_pod_residuals[index] == True:
            res_space = {}; res_mass = {}
            for ci in range(model.levelModelList[-1].coefficients.nc):
                res_space[ci] = numpy.zeros(model.levelModelList[-1].u[ci].dof.shape,'d')
                model.levelModelList[-1].getSpatialResidual(model.levelModelList[-1].u[ci].dof,res_space[ci])
                res_mass[ci] = numpy.zeros(model.levelModelList[-1].u[ci].dof.shape,'d')
                model.levelModelList[-1].getMassResidual(model.levelModelList[-1].u[ci].dof,res_space[ci])
            model.levelModelList[-1].archiveFiniteElementResiduals(self.ar[index],self.tnList[0],self.tCount,res_space,res_name_base='spatial_residual')
            model.levelModelList[-1].archiveFiniteElementResiduals(self.ar[index],self.tnList[0],self.tCount,res_mass,res_name_base='mass_residual')
    ##save model's solution values to archive
    def archiveSolution(self,model,index,t=None):
        if self.archiveFlag == ArchiveFlags.UNDEFINED:
            return
        if t == None:
            t = self.systemStepController.t_system

        log("Writing mesh header for  model = "+model.name+" at time t="+str(t),level=3)
        log("Writing solution for  model = "+model.name,level=3)
        if self.so.useOneArchive:
            if index==0:
                self.femSpaceWritten={}
            model.levelModelList[-1].archiveFiniteElementSolutions(self.ar[index],t,self.tCount,
                                                                   initialPhase=False,
                                                                   writeVectors=True,meshChanged=True,femSpaceWritten=self.femSpaceWritten,
                                                                   writeVelocityPostProcessor=self.opts.writeVPP)
        else:
            model.levelModelList[-1].archiveFiniteElementSolutions(self.ar[index],t,self.tCount,
                                                                   initialPhase=False,
                                                                   writeVectors=True,meshChanged=True,
                                                                   writeVelocityPostProcessor=self.opts.writeVPP)
        model.levelModelList[-1].archiveAnalyticalSolutions(self.ar[index],self.pList[index].analyticalSolution,
                                                            t,
                                                            self.tCount)
        #uses values in simFlags['storeQuantities']
        #q dictionary
        if self.archive_q[index] == True:
            scalarKeys = model.simTools.getScalarElementStorageKeys(model,t)
            vectorKeys = model.simTools.getVectorElementStorageKeys(model,t)
            tensorKeys = model.simTools.getTensorElementStorageKeys(model,t)
            model.levelModelList[-1].archiveElementQuadratureValues(self.ar[index],t,self.tCount,
                                                                    scalarKeys=scalarKeys,vectorKeys=vectorKeys,tensorKeys=tensorKeys,
                                                                    initialPhase=False,meshChanged=True)

        #ebq_global dictionary
        if self.archive_ebq_global[index] == True:
            scalarKeys = model.simTools.getScalarElementBoundaryStorageKeys(model,t)
            vectorKeys = model.simTools.getVectorElementBoundaryStorageKeys(model,t)
            tensorKeys = model.simTools.getTensorElementBoundaryStorageKeys(model,t)
            model.levelModelList[-1].archiveElementBoundaryQuadratureValues(self.ar[index],t,self.tCount,
                                                                            scalarKeys=scalarKeys,vectorKeys=vectorKeys,tensorKeys=tensorKeys,
                                                                            initialPhase=False,meshChanged=True)
        if self.archive_ebqe[index] == True:
            #ebqe dictionary
            scalarKeys = model.simTools.getScalarExteriorElementBoundaryStorageKeys(model,t)
            vectorKeys = model.simTools.getVectorExteriorElementBoundaryStorageKeys(model,t)
            tensorKeys = model.simTools.getTensorExteriorElementBoundaryStorageKeys(model,t)
            model.levelModelList[-1].archiveExteriorElementBoundaryQuadratureValues(self.ar[index],t,self.tCount,
                                                                                    scalarKeys=scalarKeys,vectorKeys=vectorKeys,tensorKeys=tensorKeys,
                                                                                    initialPhase=False,meshChanged=True)
        #for nonlinear POD 
        if self.archive_pod_residuals[index] == True:
            res_space = {}; res_mass = {}
            for ci in range(model.levelModelList[-1].coefficients.nc):
                res_space[ci] = numpy.zeros(model.levelModelList[-1].u[ci].dof.shape,'d')
                model.levelModelList[-1].getSpatialResidual(model.levelModelList[-1].u[ci].dof,res_space[ci])
                res_mass[ci] = numpy.zeros(model.levelModelList[-1].u[ci].dof.shape,'d')
                model.levelModelList[-1].getMassResidual(model.levelModelList[-1].u[ci].dof,res_mass[ci])
            model.levelModelList[-1].archiveFiniteElementResiduals(self.ar[index],t,self.tCount,res_space,res_name_base='spatial_residual')
            model.levelModelList[-1].archiveFiniteElementResiduals(self.ar[index],t,self.tCount,res_mass,res_name_base='mass_residual')
                
        if not self.opts.cacheArchive:
            self.ar[index].sync()
    ## clean up archive
    def closeArchive(self,model,index):
        if self.archiveFlag == None:
            return
        if self.so.useOneArchive:
            if index==0:
                log("Closing solution archive for "+self.so.name)
                self.ar[index].close()
        else:
            log("Closing solution archive for "+model.name)
            self.ar[index].close()

    def initializeViewSolution(self,model):
        """
        """
        model.viewer.preprocess(model,model.stepController.t_model_last)
        model.simTools.preprocess(model,model.stepController.t_model_last)

    ## run time visualization for modela
    def viewSolution(self,model,initialCondition=False):
        """

        """
        #mwf looking at last solns
        if (model.viewer.viewerType != 'matlab' or model.stepController.t_model_last <= self.tnList[0] or
            model.stepController.t_model_last >= self.tnList[-1]):
            model.viewer.processTimeLevel(model,model.stepController.t_model_last)
            model.simTools.processTimeLevel(model,model.stepController.t_model_last)


    ## clean up runtime visualization
    def finalizeViewSolution(self,model):
        model.viewer.postprocess(model,model.stepController.t_model_last)
        model.simTools.postprocess(model,model.stepController.t_model_last)
