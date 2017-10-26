"""
A hierarchy of classes for managing comlete numerical solution implementations

.. inheritance-diagram:: proteus.NumericalSolution
   :parts: 1
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

from .Profiling import logEvent

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

    def __init__(self,so,pList,nList,sList,opts,simFlagsList=None):
        import Comm
        comm=Comm.get()
        self.comm=comm
        message = "Initializing NumericalSolution for "+so.name+"\n System includes: \n"
        for p in pList:
            message += p.name+"\n"
        logEvent(message)
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

        logEvent("Setting Archiver(s)")

        if so.useOneArchive:
            self.femSpaceWritten={}
            tmp  = Archiver.XdmfArchive(opts.dataDir,so.name,useTextArchive=opts.useTextArchive,
                                        gatherAtClose=opts.gatherArchive,hotStart=opts.hotStart,
                                        useGlobalXMF=(not opts.subdomainArchives),
                                        global_sync=opts.global_sync)
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
        if simFlagsList is not None:
            assert len(simFlagsList) == len(self.pList), "len(simFlagsList) = %s should be %s " % (len(simFlagsList),len(self.pList))
            for index in range(len(self.pList)):
                if simFlagsList[index].has_key('storeQuantities'):
                    for quant in filter(lambda a: a is not None,simFlagsList[index]['storeQuantities']):
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
                            logEvent("Warning Numerical Solution storeQuantity = %s not recognized won't archive" % quant)
                    #
                #
            #
        #
        logEvent("Setting up MultilevelMesh")
        mlMesh_nList = []
        if so.useOneMesh:
            logEvent("Building one multilevel mesh for all models")
            nListForMeshGeneration=[nList[0]]
            pListForMeshGeneration=[pList[0]]
        else:
            logEvent("Building seperate meshes for each model")
            nListForMeshGeneration=nList
            pListForMeshGeneration=pList

        for p,n in zip(pListForMeshGeneration,nListForMeshGeneration):
            if opts.hotStart:
                p.genMesh = False
                logEvent("Hotstarting, using existing mesh "+p.name)
            else:
                logEvent("Generating mesh for "+p.name)
            #support for old-style domain input
            if p.domain is None:
                if p.nd == 1:
                    p.domain = Domain.RectangularDomain(L=p.L[:1],
                                                        x=p.x0[:1],
                                                        name=p.name)
                elif p.nd == 2:
                    if p.polyfile is not None:
                        p.domain = Domain.PlanarStraightLineGraphDomain(fileprefix=p.polyfile,name=p.polyfile)
                    elif p.meshfile != None:
                        p.domain = Domain.Mesh2DMDomain(p.meshfile)
                    else:
                        p.domain = Domain.RectangularDomain(L=p.L[:2],
                                                            x=p.x0[:2],
                                                            name=p.name)
                elif p.nd == 3:
                    if p.polyfile is not None:
                        p.domain = Domain.PiecewiseLinearComplexDomain(fileprefix=p.polyfile,name=p.polyfile)
                    elif p.meshfile is not None:
                        p.domain = Domain.Mesh3DMDomain(p.meshfile)
                    else:
                        p.domain = Domain.RectangularDomain(L=p.L[:3],
                                                            x=p.x0[:3],
                                                            name=p.name)
                else:
                    raise RuntimeError("No support for domains in more than three dimensions")
            #now generate meshes, could move to Domain and use polymorphism or MeshTools
            if isinstance(p.domain,Domain.RectangularDomain):
                if p.domain.nd == 1:
                    mlMesh = MeshTools.MultilevelEdgeMesh(n.nn, 1, 1,
                                                          p.domain.x[0], 0.0, 0.0,
                                                          p.domain.L[0], 1.0, 1.0,
                                                          refinementLevels=n.nLevels,
                                                          nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                          parallelPartitioningType=n.parallelPartitioningType)
                elif p.domain.nd == 2:
                    if (n.nnx == n.nny is None):
                        nnx = nny = n.nn
                    else:
                        nnx = n.nnx
                        nny = n.nny
                    logEvent("Building %i x %i rectangular mesh for %s" % (nnx,nny,p.name))

                    if not hasattr(n,'quad'):
                        n.quad = False

                    if (n.quad):
                        mlMesh = MeshTools.MultilevelQuadrilateralMesh(nnx,nny,1,
                                                                       p.domain.x[0], p.domain.x[1], 0.0,
                                                                       p.domain.L[0],p.domain.L[1],1,
                                                                       refinementLevels=n.nLevels,
                                                                       nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                                       parallelPartitioningType=n.parallelPartitioningType)
                    else:
                        mlMesh = MeshTools.MultilevelTriangularMesh(nnx,nny,1,
                                                                    p.domain.x[0], p.domain.x[1], 0.0,
                                                                    p.domain.L[0],p.domain.L[1],1,
                                                                    refinementLevels=n.nLevels,
                                                                    nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                                    parallelPartitioningType=n.parallelPartitioningType)

                elif p.domain.nd == 3:
                    if (n.nnx == n.nny == n.nnz  is None):
                        nnx = nny = nnz = n.nn
                    else:
                        nnx = n.nnx
                        nny = n.nny
                        nnz = n.nnz
                    logEvent("Building %i x %i x %i rectangular mesh for %s" % (nnx,nny,nnz,p.name))

                    if not hasattr(n,'hex'):
                        n.hex = False

                    if not hasattr(n,'NURBS'):
                        n.NURBS = False

                    if (n.NURBS):
                        mlMesh = MeshTools.MultilevelNURBSMesh(nnx,nny,nnz,
                                                               n.px,n.py,n.pz,
                                                               p.domain.x[0], p.domain.x[1], p.domain.x[2],
                                                               p.domain.L[0], p.domain.L[1], p.domain.L[2],
                                                               refinementLevels=n.nLevels,
                                                               nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                               parallelPartitioningType=n.parallelPartitioningType)
                    elif (n.hex):
                        if not hasattr(n,'px'):
                            n.px=0
                            n.py=0
                            n.pz=0
                        mlMesh = MeshTools.MultilevelHexahedralMesh(nnx, nny, nnz,
                                                                    n.px,n.py,n.pz,
                                                                    p.domain.x[0], p.domain.x[1], p.domain.x[2],
                                                                    p.domain.L[0], p.domain.L[1], p.domain.L[2],
                                                                    refinementLevels=n.nLevels,
                                                                    nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                                    parallelPartitioningType=n.parallelPartitioningType)
                    else :
                        mlMesh = MeshTools.MultilevelTetrahedralMesh(nnx, nny, nnz,
                                                                     p.domain.x[0], p.domain.x[1], p.domain.x[2],
                                                                     p.L[0], p.L[1], p.L[2],
                                                                     refinementLevels=n.nLevels,
                                                                     nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                                     parallelPartitioningType=n.parallelPartitioningType)

            elif isinstance(p.domain,Domain.PlanarStraightLineGraphDomain):
                if p.domain.use_gmsh is True:
                    if comm.isMaster() and (p.genMesh or not (os.path.exists(p.domain.geofile+".ele") and
                                                              os.path.exists(p.domain.geofile+".node") and
                                                              os.path.exists(p.domain.geofile+".edge"))):
                        logEvent("Running gmsh to generate 2D mesh for "+p.name,level=1)
                        gmsh_cmd = "time gmsh {0:s} -v 10 -2 -o {1:s} -format msh".format(p.domain.geofile+".geo", p.domain.geofile+".msh")
                        logEvent("Calling gmsh on rank 0 with command %s" % (gmsh_cmd,))
                        check_call(gmsh_cmd, shell=True)
                        logEvent("Done running gmsh; converting to triangle")
                        MeshTools.msh2simplex(fileprefix=p.domain.geofile, nd=2)

                    comm.barrier()
                    mesh = MeshTools.TriangularMesh()
                    mlMesh = MeshTools.MultilevelTriangularMesh(0,0,0,skipInit=True,
                                                                nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                                parallelPartitioningType=n.parallelPartitioningType)
                    logEvent("NAHeader GridRefinements %i" % (n.nLevels) )
                    logEvent("Generating %i-level mesh from coarse Triangle mesh" % (n.nLevels,))
                    logEvent("Generating coarse global mesh from Triangle files")
                    mesh.generateFromTriangleFiles(filebase=p.domain.geofile,base=1)
                    logEvent("Generating partitioned %i-level mesh from coarse global Triangle mesh" % (n.nLevels,))
                    mlMesh.generateFromExistingCoarseMesh(mesh,n.nLevels,
                                                          nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                          parallelPartitioningType=n.parallelPartitioningType)
                else:
                    logEvent("Calling Triangle to generate 2D mesh for"+p.name)
                    tmesh = TriangleTools.TriangleBaseMesh(baseFlags=n.triangleOptions,
                                                           nbase=1,
                                                           verbose=10)
                    if comm.isMaster() and p.genMesh:
                        tmesh.readFromPolyFile(p.domain.polyfile)
                        tmesh.writeToFile(p.domain.polyfile)
                    comm.barrier()

                    mesh = MeshTools.TriangularMesh()
                    mesh.generateFromTriangleFiles(filebase=p.domain.polyfile,base=1)

                    mlMesh = MeshTools.MultilevelTriangularMesh(0,0,0,skipInit=True,
                                                                nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                                parallelPartitioningType=n.parallelPartitioningType)
                    logEvent("Generating %i-level mesh from coarse Triangle mesh" % (n.nLevels,))
                    mlMesh.generateFromExistingCoarseMesh(mesh,n.nLevels,
                                                        nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                        parallelPartitioningType=n.parallelPartitioningType)

            elif isinstance(p.domain,Domain.PiecewiseLinearComplexDomain):
                from subprocess import call
                import sys
                if p.domain.use_gmsh is True:
                    fileprefix = p.domain.geofile
                else:
                    fileprefix = p.domain.polyfile
                if comm.rank() == 0 and (p.genMesh or not (os.path.exists(fileprefix+".ele") and
                                                           os.path.exists(fileprefix+".node") and
                                                           os.path.exists(fileprefix+".face"))):
                    if p.domain.use_gmsh is True:
                        logEvent("Running gmsh to generate 3D mesh for "+p.name,level=1)
                        gmsh_cmd = "time gmsh {0:s} -v 10 -3 -o {1:s} -format msh".format(fileprefix+'.geo', p.domain.geofile+'.msh')
                        logEvent("Calling gmsh on rank 0 with command %s" % (gmsh_cmd,))
                        check_call(gmsh_cmd, shell=True)
                        logEvent("Done running gmsh; converting to tetgen")
                        MeshTools.msh2simplex(fileprefix=fileprefix, nd=3)
                        check_call("tetgen -Vfeen {0:s}.ele".format(fileprefix), shell=True)
                    else:
                        logEvent("Running tetgen to generate 3D mesh for "+p.name, level=1)
                        tetcmd = "tetgen -{0} {1}.poly".format(n.triangleOptions, fileprefix)
                        logEvent("Calling tetgen on rank 0 with command %s" % (tetcmd,))
                        check_call(tetcmd, shell=True)
                        logEvent("Done running tetgen")
                    check_call("mv {0:s}.1.ele {0:s}.ele".format(fileprefix), shell=True)
                    check_call("mv {0:s}.1.node {0:s}.node".format(fileprefix), shell=True)
                    check_call("mv {0:s}.1.face {0:s}.face".format(fileprefix), shell=True)
                    try:
                        check_call("mv {0:s}.1.neigh {0:s}.neigh".format(fileprefix), shell=True)
                    except:
                        logEvent("Warning: couldn't move {0:s}.1.neigh".format(fileprefix))
                        pass
                    try:
                        logEvent("Warning: couldn't move {0:s}.1.edge".format(fileprefix))
                        check_call("mv {0:s}.1.edge {0:s}.edge".format(fileprefix), shell=True)
                    except:
                        pass
                comm.barrier()
                logEvent("Initializing mesh and MultilevelMesh")
                nbase = 1
                mesh=MeshTools.TetrahedralMesh()
                mlMesh = MeshTools.MultilevelTetrahedralMesh(0,0,0,skipInit=True,
                                                             nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                             parallelPartitioningType=n.parallelPartitioningType)
                if opts.generatePartitionedMeshFromFiles:
                    logEvent("Generating partitioned mesh from Tetgen files")
                    mlMesh.generatePartitionedMeshFromTetgenFiles(fileprefix,nbase,mesh,n.nLevels,
                                                                  nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                                  parallelPartitioningType=n.parallelPartitioningType)
                else:
                    logEvent("Generating coarse global mesh from Tetgen files")
                    mesh.generateFromTetgenFiles(fileprefix,nbase,parallel = comm.size() > 1)
                    logEvent("Generating partitioned %i-level mesh from coarse global Tetgen mesh" % (n.nLevels,))
                    mlMesh.generateFromExistingCoarseMesh(mesh,n.nLevels,
                                                          nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                          parallelPartitioningType=n.parallelPartitioningType)
            elif isinstance(p.domain,Domain.PUMIDomain):
                #ibaned: PUMI conversion #1
                if p.domain.nd == 3:
                  mesh = MeshTools.TetrahedralMesh()
                else:
                  mesh = MeshTools.TriangularMesh()
                logEvent("Converting PUMI mesh to Proteus")
                mesh.convertFromPUMI(p.domain.PUMIMesh, p.domain.faceList,
                    p.domain.regList,
                    parallel = comm.size() > 1, dim = p.domain.nd)
                if p.domain.nd == 3:
                  mlMesh = MeshTools.MultilevelTetrahedralMesh(
                      0,0,0,skipInit=True,
                      nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                      parallelPartitioningType=n.parallelPartitioningType)
                if p.domain.nd == 2:
                  mlMesh = MeshTools.MultilevelTriangularMesh(
                      0,0,0,skipInit=True,
                      nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                      parallelPartitioningType=n.parallelPartitioningType)
                logEvent("Generating %i-level mesh from PUMI mesh" % (n.nLevels,))
                if comm.size()==1:
                  mlMesh.generateFromExistingCoarseMesh(
                      mesh,n.nLevels,
                      nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                      parallelPartitioningType=n.parallelPartitioningType)
                else:
                  mlMesh.generatePartitionedMeshFromPUMI(
                      mesh,n.nLevels,
                      nLayersOfOverlap=n.nLayersOfOverlapForParallel)
            elif isinstance(p.domain,Domain.MeshTetgenDomain):
                nbase = 1
                mesh=MeshTools.TetrahedralMesh()
                logEvent("Reading coarse mesh from tetgen file")
                mlMesh = MeshTools.MultilevelTetrahedralMesh(0,0,0,skipInit=True,
                                                             nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                             parallelPartitioningType=n.parallelPartitioningType)
                if opts.generatePartitionedMeshFromFiles:
                    logEvent("Generating partitioned mesh from Tetgen files")
                    mlMesh.generatePartitionedMeshFromTetgenFiles(p.domain.meshfile,nbase,mesh,n.nLevels,
                                                                  nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                                  parallelPartitioningType=n.parallelPartitioningType)
                else:
                    logEvent("Generating coarse global mesh from Tetgen files")
                    mesh.generateFromTetgenFiles(p.domain.polyfile,nbase,parallel = comm.size() > 1)
                    logEvent("Generating partitioned %i-level mesh from coarse global Tetgen mesh" % (n.nLevels,))
                    mlMesh.generateFromExistingCoarseMesh(mesh,n.nLevels,
                                                          nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                          parallelPartitioningType=n.parallelPartitioningType)
            elif isinstance(p.domain,Domain.Mesh3DMDomain):
                mesh=MeshTools.TetrahedralMesh()
                logEvent("Reading coarse mesh from 3DM file")
                mesh.generateFrom3DMFile(p.domain.meshfile)
                mlMesh = MeshTools.MultilevelTetrahedralMesh(0,0,0,skipInit=True,
                                                             nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                             parallelPartitioningType=n.parallelPartitioningType)
                logEvent("Generating %i-level mesh from coarse 3DM mesh" % (n.nLevels,))
                mlMesh.generateFromExistingCoarseMesh(mesh,n.nLevels,
                                                      nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                      parallelPartitioningType=n.parallelPartitioningType)
            elif isinstance(p.domain,Domain.Mesh2DMDomain):
                mesh=MeshTools.TriangularMesh()
                logEvent("Reading coarse mesh from 2DM file")
                mesh.generateFrom2DMFile(p.domain.meshfile)
                mlMesh = MeshTools.MultilevelTriangularMesh(0,0,0,skipInit=True,
                                                             nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                             parallelPartitioningType=n.parallelPartitioningType)
                logEvent("Generating %i-level mesh from coarse 2DM mesh" % (n.nLevels,))
                mlMesh.generateFromExistingCoarseMesh(mesh,n.nLevels,
                                                      nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                      parallelPartitioningType=n.parallelPartitioningType)
            elif isinstance(p.domain,Domain.MeshHexDomain):
                mesh=MeshTools.HexahedralMesh()
                logEvent("Reading coarse mesh from file")
                mesh.generateFromHexFile(p.domain.meshfile)
                mlMesh = MeshTools.MultilevelHexahedralMesh(0,0,0,skipInit=True,
                                                             nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                             parallelPartitioningType=n.parallelPartitioningType)
                logEvent("Generating %i-level mesh from coarse mesh" % (n.nLevels,))
                mlMesh.generateFromExistingCoarseMesh(mesh,n.nLevels,
                                                      nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                      parallelPartitioningType=n.parallelPartitioningType)
            elif isinstance(p.domain,Domain.GMSH_3D_Domain):
                from subprocess import call
                import sys
                if comm.rank() == 0 and (p.genMesh or not (os.path.exists(p.domain.polyfile+".ele") and
                                                           os.path.exists(p.domain.polyfile+".node") and
                                                           os.path.exists(p.domain.polyfile+".face"))):
                    logEvent("Running gmsh to generate 3D mesh for "+p.name,level=1)
                    gmsh_cmd = "time gmsh {0:s} -v 10 -3 -o {1:s}  -format mesh  -clmax {2:f} -clscale {2:f}".format(p.domain.geofile, p.domain.name+".mesh", p.domain.he)

                    logEvent("Calling gmsh on rank 0 with command %s" % (gmsh_cmd,))

                    check_call(gmsh_cmd, shell=True)

                    logEvent("Done running gmsh; converting to tetgen")

                    gmsh2tetgen_cmd = "gmsh2tetgen {0} {1:f} {2:d} {3:d} {4:d}".format(
                        p.domain.name+".mesh",
                        p.domain.length_scale,
                        p.domain.permute_dims[0]+1,
                        p.domain.permute_dims[1]+1,
                        p.domain.permute_dims[2]+1)

                    check_call(gmsh2tetgen_cmd, shell=True)
                    check_call("tetgen -Vfeen %s.ele" % ("mesh",), shell=True)
                    check_call("mv %s.1.ele %s.ele" % ("mesh","mesh"), shell=True)
                    check_call("mv %s.1.node %s.node" % ("mesh","mesh"), shell=True)
                    check_call("mv %s.1.face %s.face" % ("mesh","mesh"), shell=True)
                    check_call("mv %s.1.neigh %s.neigh" % ("mesh","mesh"), shell=True)
                    check_call("mv %s.1.edge %s.edge" % ("mesh","mesh"), shell=True)
                    elefile  = "mesh.ele"
                    nodefile = "mesh.node"
                    facefile = "mesh.face"
                    edgefile = "mesh.edge"
                    assert os.path.exists(elefile), "no mesh.ele"
                    tmp = "%s.ele" % p.domain.polyfile
                    os.rename(elefile,tmp)
                    assert os.path.exists(tmp), "no .ele"
                    assert os.path.exists(nodefile), "no mesh.node"
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
                logEvent("Initializing mesh and MultilevelMesh")
                nbase = 1
                mesh=MeshTools.TetrahedralMesh()
                mlMesh = MeshTools.MultilevelTetrahedralMesh(0,0,0,skipInit=True,
                                                             nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                             parallelPartitioningType=n.parallelPartitioningType)
                if opts.generatePartitionedMeshFromFiles:
                    logEvent("Generating partitioned mesh from Tetgen files")
                    mlMesh.generatePartitionedMeshFromTetgenFiles(p.domain.polyfile,nbase,mesh,n.nLevels,
                                                                  nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                                  parallelPartitioningType=n.parallelPartitioningType)
                else:
                    logEvent("Generating coarse global mesh from Tetgen files")
                    mesh.generateFromTetgenFiles(p.domain.polyfile,nbase,parallel = comm.size() > 1)
                    logEvent("Generating partitioned %i-level mesh from coarse global Tetgen mesh" % (n.nLevels,))
                    mlMesh.generateFromExistingCoarseMesh(mesh,n.nLevels,
                                                          nLayersOfOverlap=n.nLayersOfOverlapForParallel,
                                                          parallelPartitioningType=n.parallelPartitioningType)
            mlMesh_nList.append(mlMesh)
            if opts.viewMesh:
                logEvent("Attempting to visualize mesh")
                try:
                    from proteusGraphical import vtkViewers
                    vtkViewers.ViewMesh(mlMesh.meshList[0],viewMaterialTypes=True)
                    vtkViewers.ViewBoundaryMesh(mlMesh.meshList[0],viewBoundaryMaterialTypes=True)
                except:
                    logEvent("NumericalSolution ViewMesh failed for coarse mesh")
            for l in range(n.nLevels):
                try:
                    logEvent(mlMesh.meshList[l].meshInfo())
                except:
                    logEvent("meshInfo() method not implemented for this mesh type")
                if opts.viewMesh and opts.viewLevels and l > 0:
                    logEvent("Attempting to visualize mesh")
                    try:
                        from proteusGraphical import vtkViewers
                        vtkViewers.ViewMesh(mlMesh.meshList[l],title="mesh level %s " % l,
                                            viewMaterialTypes=True)
                        vtkViewers.ViewBoundaryMesh(mlMesh.meshList[l],title="boundary mesh level %s " % l,
                                                    viewBoundaryMaterialTypes=True)
                    except:
                        logEvent("NumericalSolution ViewMesh failed for mesh level %s" % l)
        if so.useOneMesh:
            for p in pList[1:]: mlMesh_nList.append(mlMesh)
            try:
                if (nList[0].MeshAdaptMesh.size_field_config() == 'isotropicProteus'):
                    mlMesh.meshList[0].subdomainMesh.size_field = numpy.ones((mlMesh.meshList[0].subdomainMesh.nNodes_global,1),'d')*1.0e-1
                if (nList[0].MeshAdaptMesh.size_field_config() == 'anisotropicProteus'):
                    mlMesh.meshList[0].subdomainMesh.size_scale = numpy.ones((mlMesh.meshList[0].subdomainMesh.nNodes_global,3),'d')
                    mlMesh.meshList[0].subdomainMesh.size_frame = numpy.ones((mlMesh.meshList[0].subdomainMesh.nNodes_global,9),'d')
            except:
                pass
        Profiling.memory("Mesh")
        from collections import OrderedDict
        self.modelSpinUp = OrderedDict()
        for p in pList:
            p.coefficients.opts = self.opts
            if p.coefficients.sdInfo == {}:
                for ci,ckDict in p.coefficients.diffusion.iteritems():
                    for ck in ckDict.keys():
                        if not p.coefficients.sdInfo.has_key((ci,ck)):
                            p.coefficients.sdInfo[(ci,ck)] = (numpy.arange(start=0,stop=p.nd**2+1,step=p.nd,dtype='i'),
                                                              numpy.array([range(p.nd) for row in range(p.nd)],dtype='i').flatten())
                            logEvent("Numerical Solution Sparse diffusion information key "+`(ci,ck)`+' = '+`p.coefficients.sdInfo[(ci,ck)]`)
        self.sList = sList
        self.mlMesh_nList = mlMesh_nList
        self.allocateModels()
        #collect models to be used for spin up
        for index in so.modelSpinUpList:
            self.modelSpinUp[index] = self.modelList[index]
        logEvent("Finished setting up models and solvers")
        if self.opts.save_dof:
            for m in self.modelList:
                for lm in m.levelModelList:
                    for ci in range(lm.coefficients.nc):
                        lm.u[ci].dof_last = lm.u[ci].dof.copy()
        self.archiveFlag= so.archiveFlag
        logEvent("Setting up SimTools for "+p.name)
        self.simOutputList = []
        self.auxiliaryVariables = {}
        if self.simFlagsList is not None:
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
        logEvent(Profiling.memory("NumericalSolution memory",className='NumericalSolution',memSaved=memBase))
        if so.tnList is None:
            logEvent("Building tnList from model = "+pList[0].name+" nDTout = "+`nList[0].nDTout`)
            self.tnList=[float(n)*nList[0].T/float(nList[0].nDTout)
                         for n in range(nList[0].nDTout+1)]
        else:
            logEvent("Using tnList from so = "+so.name)
            self.tnList = so.tnList
        logEvent("Time sequence"+`self.tnList`)
        logEvent("NAHeader Num Time Steps "+`len(self.tnList)-1`)
        logEvent("Setting "+so.name+" systemStepController to object of type "+str(so.systemStepControllerType))
        self.systemStepController = so.systemStepControllerType(self.modelList,stepExact=so.systemStepExact)
        self.systemStepController.setFromOptions(so)
        logEvent("Finished NumericalSolution initialization")

    def allocateModels(self):
        self.modelList=[]
        self.lsList=[]
        self.nlsList=[]
        
        for p,n,s,mlMesh,index \
            in zip(self.pList,self.nList,self.sList,self.mlMesh_nList,range(len(self.pList))):

            if self.so.needEBQ_GLOBAL:
                n.needEBQ_GLOBAL = True

            if self.so.needEBQ:
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

            logEvent("Setting up MultilevelTransport for "+p.name)

            model \
                = Transport.MultilevelTransport(p,
                                                n,
                                                mlMesh,
                                                OneLevelTransportType=p.LevelModelType)
            self.modelList.append(model)

            model.name = p.name
            logEvent("Setting "+model.name+" stepController to "+str(n.stepController))
            model.stepController = n.stepController(model,n)

            Profiling.memory("MultilevelTransport for "+p.name)
            logEvent("Setting up MultilevelLinearSolver for"+p.name)

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
                solver_options_prefix=linear_solver_options_prefix,
                computeEigenvalues = n.computeEigenvalues,
                linearSmootherOptions = n.linearSmootherOptions)
            self.lsList.append(multilevelLinearSolver)
            Profiling.memory("MultilevelLinearSolver for "+p.name)
            logEvent("Setting up MultilevelNonLinearSolver for "+p.name)
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

    def PUMI2Proteus(self,mesh):
        p0 = self.pList[0] #This can probably be cleaned up somehow
        n0 = self.nList[0]
        logEvent("Generating %i-level mesh from PUMI mesh" % (n0.nLevels,))
        if p0.domain.nd == 3:
          mlMesh = MeshTools.MultilevelTetrahedralMesh(
              0,0,0,skipInit=True,
              nLayersOfOverlap=n0.nLayersOfOverlapForParallel,
              parallelPartitioningType=n0.parallelPartitioningType)
        if p0.domain.nd == 2:
          mlMesh = MeshTools.MultilevelTriangularMesh(
              0,0,0,skipInit=True,
              nLayersOfOverlap=n0.nLayersOfOverlapForParallel,
              parallelPartitioningType=n0.parallelPartitioningType)
        if self.comm.size()==1:
            mlMesh.generateFromExistingCoarseMesh(
                mesh,n0.nLevels,
                nLayersOfOverlap=n0.nLayersOfOverlapForParallel,
                parallelPartitioningType=n0.parallelPartitioningType)
        else:
            mlMesh.generatePartitionedMeshFromPUMI(
                mesh,n0.nLevels,
                nLayersOfOverlap=n0.nLayersOfOverlapForParallel)
        self.mlMesh_nList=[]
        for p in self.pList:
            self.mlMesh_nList.append(mlMesh)
        if (p0.domain.PUMIMesh.size_field_config() == "isotropicProteus"):
            mlMesh.meshList[0].subdomainMesh.size_field = numpy.ones((mlMesh.meshList[0].subdomainMesh.nNodes_global,1),'d')*1.0e-1
        if (p0.domain.PUMIMesh.size_field_config() == 'anisotropicProteus'):
            mlMesh.meshList[0].subdomainMesh.size_scale = numpy.ones((mlMesh.meshList[0].subdomainMesh.nNodes_global,3),'d')
            mlMesh.meshList[0].subdomainMesh.size_frame = numpy.ones((mlMesh.meshList[0].subdomainMesh.nNodes_global,9),'d')

        #may want to trigger garbage collection here
        modelListOld = self.modelList
        logEvent("Allocating models on new mesh")
        self.allocateModels()
        logEvent("Attach auxiliary variables to new models")
        #(cut and pasted from init, need to cleanup)
        self.simOutputList = []
        self.auxiliaryVariables = {}
        self.newAuxiliaryVariables = {}
        if self.simFlagsList is not None:
            for p, n, simFlags, model, index in zip(
                    self.pList,
                    self.nList,
                    self.simFlagsList,
                    self.modelList,
                    range(len(self.pList))):
                self.simOutputList.append(
                    SimTools.SimulationProcessor(
                        flags=simFlags,
                        nLevels=n.nLevels,
                        pFile=p,
                        nFile=n,
                        analyticalSolution=p.analyticalSolution))
                model.simTools = self.simOutputList[-1]
                
                #Code to refresh attached gauges. The goal is to first purge 
                #existing point gauge node associations as that may have changed
                #If there is a line gauge, then all the points must be deleted
                #and remade.
                from collections import OrderedDict
                for av in n.auxiliaryVariables:
                  if hasattr(av,'adapted'):
                    av.adapted=True
                    for point, l_d in av.points.iteritems():
                      if 'nearest_node' in l_d:
                        l_d.pop('nearest_node')
                    if(av.isLineGauge or av.isLineIntegralGauge): #if line gauges, need to remove all points
                      av.points = OrderedDict()
                    if(av.isGaugeOwner):
                      if(self.comm.rank()==0 and not av.file.closed):
                        av.file.close()
                      for item in av.pointGaugeVecs:
                        item.destroy()
                      for item in av.pointGaugeMats:
                        item.destroy()
                      for item in av.dofsVecs:
                        item.destroy()

                      av.pointGaugeVecs = []
                      av.pointGaugeMats = []
                      av.dofsVecs = []
                      av.field_ids=[]
                      av.isGaugeOwner=False
                ##reinitialize auxiliaryVariables
                self.auxiliaryVariables[model.name]= [av.attachModel(model,self.ar[index]) for av in n.auxiliaryVariables]
        else:
            for p,n,s,model,index in zip(
                    self.pList,
                    self.nList,
                    self.sList,
                    self.modelList,
                    range(len(self.pList))):
                self.simOutputList.append(SimTools.SimulationProcessor(pFile=p,nFile=n))
                model.simTools = self.simOutputList[-1]
                model.viewer = Viewers.V_base(p,n,s)
                self.auxiliaryVariables[model.name]= [av.attachModel(model,self.ar[index]) for av in n.auxiliaryVariables]
        for avList in self.auxiliaryVariables.values():
            for av in avList:
                av.attachAuxiliaryVariables(self.auxiliaryVariables)
        logEvent("Transfering fields from PUMI to Proteus")
        for m in self.modelList:
          for lm in m.levelModelList:
            coef = lm.coefficients
            if coef.vectorComponents is not None:
              vector=numpy.zeros((lm.mesh.nNodes_global,3),'d')
              p0.domain.PUMIMesh.transferFieldToProteus(
                     coef.vectorName, vector)
              for vci in range(len(coef.vectorComponents)):
                lm.u[coef.vectorComponents[vci]].dof[:] = vector[:,vci]
              del vector
            for ci in range(coef.nc):
              if coef.vectorComponents is None or \
                 ci not in coef.vectorComponents:
                scalar=numpy.zeros((lm.mesh.nNodes_global,1),'d')
                p0.domain.PUMIMesh.transferFieldToProteus(
                    coef.variableNames[ci], scalar)
                lm.u[ci].dof[:] = scalar[:,0]
                del scalar
        logEvent("Attaching models on new mesh to each other")
        for m,ptmp,mOld in zip(self.modelList, self.pList, modelListOld):
            for lm, lu, lr, lmOld in zip(m.levelModelList, m.uList, m.rList,mOld.levelModelList):
                save_dof=[]
                for ci in range(lm.coefficients.nc):
                    save_dof.append( lm.u[ci].dof.copy())
                    lm.u[ci].dof_last = lm.u[ci].dof.copy()
                lm.setFreeDOF(lu)
                for ci in range(lm.coefficients.nc):
                    assert((save_dof[ci] == lm.u[ci].dof).all())
                lm.calculateSolutionAtQuadrature()
                lm.timeIntegration.tLast = lmOld.timeIntegration.tLast
                lm.timeIntegration.t = lmOld.timeIntegration.t
                lm.timeIntegration.dt = lmOld.timeIntegration.dt
                assert(lmOld.timeIntegration.tLast == lm.timeIntegration.tLast)
                assert(lmOld.timeIntegration.t == lm.timeIntegration.t)
                assert(lmOld.timeIntegration.dt == lm.timeIntegration.dt)
            m.stepController.dt_model = mOld.stepController.dt_model
            m.stepController.t_model = mOld.stepController.t_model
            m.stepController.t_model_last = mOld.stepController.t_model_last
            m.stepController.substeps = mOld.stepController.substeps
        # logEvent("Evaluating residuals and time integration")
        for m,ptmp,mOld in zip(self.modelList, self.pList, modelListOld):
            logEvent("Attaching models to model "+ptmp.name)
            m.attachModels(self.modelList)
        logEvent("Evaluating residuals and time integration")
        for m,ptmp,mOld in zip(self.modelList, self.pList, modelListOld):
            for lm, lu, lr, lmOld in zip(m.levelModelList, m.uList, m.rList, mOld.levelModelList):
                lm.timeTerm=True
                lm.getResidual(lu,lr)
                lm.timeIntegration.initializeTimeHistory(resetFromDOF=True)
                lm.initializeTimeHistory()
                lm.timeIntegration.initializeSpaceHistory()
                lm.getResidual(lu,lr)
                # assert(lmOld.timeIntegration.tLast == lm.timeIntegration.tLast)
                # assert(lmOld.timeIntegration.t == lm.timeIntegration.t)
                # assert(lmOld.timeIntegration.dt == lm.timeIntegration.dt)
                #lm.coefficients.evaluate(self.t_stepSequence,lm.q)
                #lm.coefficients.evaluate(self.t_stepSequence,lm.ebqe)
                #lm.timeIntegration.calculateElementCoefficients(lm.q)
            assert(m.stepController.dt_model == mOld.stepController.dt_model)
            assert(m.stepController.t_model == mOld.stepController.t_model)
            assert(m.stepController.t_model_last == mOld.stepController.t_model_last)
            logEvent("Initializing time history for model step controller")
            m.stepController.initializeTimeHistory()
        p0.domain.initFlag=True #For next step to take initial conditions from solution, only used on restarts
        self.systemStepController.modelList = self.modelList
        self.systemStepController.exitModelStep = {}
        self.systemStepController.controllerList = []
        for model in self.modelList:
            self.systemStepController.exitModelStep[model] = False
            if model.levelModelList[-1].timeIntegration.isAdaptive:
                self.systemStepController.controllerList.append(model)
                self.systemStepController.maxFailures = model.stepController.maxSolverFailures
        self.systemStepController.choose_dt_system()
        # for m,ptmp,mOld in zip(self.modelList, self.pList, modelListOld):
        #     for lm, lu, lr, lmOld in zip(m.levelModelList, m.uList, m.rList, mOld.levelModelList):
        #         assert(lmOld.timeIntegration.tLast == lm.timeIntegration.tLast)
        #         assert(lmOld.timeIntegration.t == lm.timeIntegration.t)
        #         assert(lmOld.timeIntegration.dt == lm.timeIntegration.dt)
        #     assert(m.stepController.dt_model == mOld.stepController.dt_model)
        #     assert(m.stepController.t_model == mOld.stepController.t_model)
        #     assert(m.stepController.t_model_last == mOld.stepController.t_model_last)
        if self.archiveFlag == ArchiveFlags.EVERY_SEQUENCE_STEP:
            #hack for archiving initial solution on adapted mesh
            self.tCount+=1
            for index,model in enumerate(self.modelList):
                self.archiveSolution(
                    model,
                    index,
                    self.systemStepController.t_system_last+1.0e-6)

    def PUMI_estimateError(self):
        """
        Estimate the error using the classical element residual method by
        Ainsworth and Oden and generates a corresponding error field.
        """

        p0 = self.pList[0]
        n0 = self.nList[0]
        adaptMeshNow = False
        #will need to move this to earlier when the mesh is created
        from proteus.MeshAdaptPUMI import MeshAdaptPUMI
        if not hasattr(p0.domain,'PUMIMesh') and not isinstance(p0.domain,Domain.PUMIDomain) and n0.adaptMesh:
            import sys
            if(self.comm.size()>1 and p0.domain.MeshOptions.parallelPartitioningType!=MeshTools.MeshParallelPartitioningTypes.element):
                sys.exit("The mesh must be partitioned by elements and NOT nodes for adaptivity functionality. Do this with: `domain.MeshOptions.setParallelPartitioningType('element')'.")
            p0.domain.PUMIMesh=n0.MeshAdaptMesh
            p0.domain.hasModel = n0.useModel
            numModelEntities = numpy.array([len(p0.domain.vertices),len(p0.domain.segments),len(p0.domain.facets),len(p0.domain.regions)]).astype("i")
            #force initialization of arrays to enable passage through to C++ code
            mesh2Model_v= numpy.asarray([[0,0]]).astype("i")
            mesh2Model_e=numpy.asarray([[0,0]]).astype("i")
            mesh2Model_b=numpy.asarray([[0,0]]).astype("i")

            segmentList = numpy.asarray([[0,0]]).astype("i")
            newFacetList = numpy.asarray([[0,0]]).astype("i")
            #only appropriate for 2D use at the moment
            if p0.domain.vertices and p0.domain.hasModel and p0.domain.nd==2:
              p0.domain.getMesh2ModelClassification(self.modelList[0].levelModelList[0].mesh)
              segmentList = numpy.asarray(p0.domain.segments).astype("i")
              #force initialize the unused arrays for proper cythonization
              import copy
              newFacetList = []
              if(not p0.domain.facets):
                p0.domain.facets = [(-1,-1)]
                newFacetList = copy.deepcopy(p0.domain.facets)
              else:
                facetList = []
                maxFacetLength = 0
                numHoles = len(p0.domain.holes)
                if(numHoles): #if there are holes, there can be multiple lists of facets
                  for i in range(numHoles,len(p0.domain.facets)):
                    for j in range(len(p0.domain.facets[i])):
                      maxFacetLength = max(maxFacetLength,len(p0.domain.facets[i][j]))
                  for i in range(numHoles,len(p0.domain.facets)):
                    facetList.append(list(p0.domain.facets[i][0]))
                    if(len(p0.domain.facets[i][0])<maxFacetLength):
                      initLength = len(p0.domain.facets[i][0])
                      lenDiff = maxFacetLength-initLength
                      for j in range(lenDiff):
                        facetList[i-numHoles].append(-1)
                else:
                  for i in range(len(p0.domain.facets)):
                    maxFacetLength = max(maxFacetLength,len(p0.domain.facets[i]))
                  for i in range(len(p0.domain.facets)):
                    facetList.append(list(p0.domain.facets[i]))
                    if(len(p0.domain.facets[i])<maxFacetLength):
                      initLength = len(p0.domain.facets[i])
                      lenDiff = maxFacetLength-initLength
                      for j in range(lenDiff):
                        facetList[i-numHoles].append(-1)

                #substitute the vertex IDs with segment IDs
                newFacetList = copy.deepcopy(facetList)
                for i in range(len(facetList)):
                  for j in range(maxFacetLength):
                    if(j==maxFacetLength-1 or facetList[i][j+1]==-1):
                      testSegment = [facetList[i][j],facetList[i][0]] 
                    else:
                      testSegment = [facetList[i][j],facetList[i][j+1]]
                    try:
                      edgIdx = p0.domain.segments.index(testSegment)
                    except ValueError:
                      edgIdx = p0.domain.segments.index(list(reversed(testSegment)))
                    newFacetList[i][j] = edgIdx
                    if(j==maxFacetLength-1 or facetList[i][j+1]==-1):
                      break
              newFacetList = numpy.asarray(newFacetList).astype("i")
              mesh2Model_v = numpy.asarray(p0.domain.meshVertex2Model).astype("i")
              mesh2Model_e = numpy.asarray(p0.domain.meshEdge2Model).astype("i")
              mesh2Model_b = numpy.asarray(p0.domain.meshBoundary2Model).astype("i")
            
            p0.domain.PUMIMesh.transferModelInfo(numModelEntities,segmentList,newFacetList,mesh2Model_v,mesh2Model_e,mesh2Model_b)
            p0.domain.PUMIMesh.reconstructFromProteus(self.modelList[0].levelModelList[0].mesh.cmesh,self.modelList[0].levelModelList[0].mesh.globalMesh.cmesh,p0.domain.hasModel)
        if (hasattr(p0.domain, 'PUMIMesh') and
            n0.adaptMesh and
            self.so.useOneMesh and 
            self.nSolveSteps%n0.adaptMesh_nSteps==0):
            logEvent("Copying coordinates to PUMI")
            p0.domain.PUMIMesh.transferFieldToPUMI("coordinates",
                self.modelList[0].levelModelList[0].mesh.nodeArray)
            if (p0.domain.PUMIMesh.size_field_config() == "isotropicProteus"):
                p0.domain.PUMIMesh.transferFieldToPUMI("proteus_size",
                                                       self.modelList[0].levelModelList[0].mesh.size_field)
            if (p0.domain.PUMIMesh.size_field_config() == 'anisotropicProteus'):
                #Insert a function to define the size_scale/size_frame fields here.
                #For a given vertex, the i-th size_scale is roughly the desired edge length along the i-th direction specified by the size_frame
                for i in range(len(self.modelList[0].levelModelList[0].mesh.size_scale)):
                  self.modelList[0].levelModelList[0].mesh.size_scale[i,0] =  1e-1
                  self.modelList[0].levelModelList[0].mesh.size_scale[i,1] =  (self.modelList[0].levelModelList[0].mesh.nodeArray[i,1]/0.584)*1e-1
                  for j in range(3):
                    for k in range(3):
                      if(j==k):
                        self.modelList[0].levelModelList[0].mesh.size_frame[i,3*j+k] = 1.0
                      else:
                        self.modelList[0].levelModelList[0].mesh.size_frame[i,3*j+k] = 0.0
                self.modelList[0].levelModelList[0].mesh.size_scale
                p0.domain.PUMIMesh.transferFieldToPUMI("proteus_sizeScale", self.modelList[0].levelModelList[0].mesh.size_scale)
                p0.domain.PUMIMesh.transferFieldToPUMI("proteus_sizeFrame", self.modelList[0].levelModelList[0].mesh.size_frame)
            p0.domain.PUMIMesh.transferFieldToPUMI("coordinates",
                self.modelList[0].levelModelList[0].mesh.nodeArray)
            logEvent("Copying DOF and parameters to PUMI")
            for m in self.modelList:
              for lm in m.levelModelList:
                coef = lm.coefficients
                if coef.vectorComponents is not None:
                  vector=numpy.zeros((lm.mesh.nNodes_global,3),'d')
                  for vci in range(len(coef.vectorComponents)):
                    vector[:,vci] = lm.u[coef.vectorComponents[vci]].dof[:]
                  p0.domain.PUMIMesh.transferFieldToPUMI(
                         coef.vectorName, vector)
                  del vector
                for ci in range(coef.nc):
                  if coef.vectorComponents is None or \
                     ci not in coef.vectorComponents:
                    scalar=numpy.zeros((lm.mesh.nNodes_global,1),'d')
                    scalar[:,0] = lm.u[ci].dof[:]
                    p0.domain.PUMIMesh.transferFieldToPUMI(
                        coef.variableNames[ci], scalar)
                    del scalar

            scalar=numpy.zeros((lm.mesh.nNodes_global,1),'d')
            # scalar[:,0] = self.modelList[0].levelModelList[0].velocityErrorNodal           
            # p0.domain.PUMIMesh.transferFieldToPUMI(
            #     'velocityError', scalar)

            # This is hardcoded for the RANS3PF to be the 4th model
            scalar[:,0] = self.modelList[4].levelModelList[0].coefficients.phi_s
            p0.domain.PUMIMesh.transferFieldToPUMI(
                'phi_s', scalar)

            del scalar
            #Get Physical Parameters
            #Can we do this in a problem-independent  way?
            rho = numpy.array([self.pList[0].rho_0,
                               self.pList[0].rho_1])
            nu = numpy.array([self.pList[0].nu_0,
                              self.pList[0].nu_1])
            g = numpy.asarray(self.pList[0].g)
            deltaT = self.tn-self.tn_last
            p0.domain.PUMIMesh.transferPropertiesToPUMI(rho,nu,g)
            del rho, nu, g

            logEvent("Estimate Error")
            sfConfig = p0.domain.PUMIMesh.size_field_config()
            if(sfConfig=="ERM"):
              errorTotal= p0.domain.PUMIMesh.get_local_error()
              
              if(p0.domain.PUMIMesh.willAdapt()):
                adaptMeshNow=True
                logEvent("Need to Adapt")
            elif(sfConfig=='interface' ):
              adaptMeshNow=True
              logEvent("Need to Adapt")
            elif(sfConfig=='meshQuality'):
              minQual = p0.domain.PUMIMesh.getMinimumQuality()
              if(minQual):
                logEvent('The quality is %f ' % minQual)
              adaptMeshNow=True
              logEvent("Need to Adapt")
            else:
              adaptMeshNow=True
              logEvent("Need to Adapt")
        return adaptMeshNow


    def PUMI_adaptMesh(self):
        """
        Uses a computed error field to construct a size field and adapts
        the mesh using SCOREC tools (a.k.a. MeshAdapt)
        """
        ##
        ## zhang-alvin's BC communication for N-S error estimation
        ##
        #  #for idx in range (0, self.modelList[0].levelModelList[0].coefficients.nc):
        #    #if idx>0:
        #    #    diff_flux = self.modelList[0].levelModelList[0].ebqe[('diffusiveFlux_bc',idx,idx)]
        #    #else:
        #    #    diff_flux = numpy.empty([2,2]) #dummy diff flux
        #    #p.domain.PUMIMesh.transferBCtagsToProteus(
        #    #    self.modelList[0].levelModelList[0].numericalFlux.isDOFBoundary[idx],
        #    #    idx,
        #    #    self.modelList[0].levelModelList[0].numericalFlux.mesh.exteriorElementBoundariesArray,
        #    #    self.modelList[0].levelModelList[0].numericalFlux.mesh.elementBoundaryElementsArray,
        #    #    diff_flux)
        #    #p.domain.PUMIMesh.transferBCtagsToProteus(
        #    #    self.modelList[0].levelModelList[0].numericalFlux.isDiffusiveFluxBoundary[idx],
        #    #    idx,
        #    #    self.modelList[0].levelModelList[0].numericalFlux.mesh.exteriorElementBoundariesArray,
        #    #    self.modelList[0].levelModelList[0].numericalFlux.mesh.elementBoundaryElementsArray,
        #    #    diff_flux)
        p0 = self.pList[0]
        n0 = self.nList[0]
        sfConfig = p0.domain.PUMIMesh.size_field_config()
        logEvent("h-adapt mesh by calling AdaptPUMIMesh")
        p0.domain.PUMIMesh.adaptPUMIMesh()

        #code to suggest adapting until error is reduced;
        #not fully baked and can lead to infinite loops of adaptation
        #if(sfConfig=="ERM"):
        #  p0.domain.PUMIMesh.get_local_error() 
        #  while(p0.domain.PUMIMesh.willAdapt()):
        #    p0.domain.PUMIMesh.adaptPUMIMesh()
        #    p0.domain.PUMIMesh.get_local_error()
        
        logEvent("Converting PUMI mesh to Proteus")
        #ibaned: PUMI conversion #2
        #TODO: this code is nearly identical to
        #PUMI conversion #1, they should be merged
        #into a function
        if p0.domain.nd == 3:
          mesh = MeshTools.TetrahedralMesh()
        else:
          mesh = MeshTools.TriangularMesh()
    
        mesh.convertFromPUMI(p0.domain.PUMIMesh,
                             p0.domain.faceList,
                             p0.domain.regList,
                             parallel = self.comm.size() > 1,
                             dim = p0.domain.nd)
        self.PUMI2Proteus(mesh)
      ##chitak end Adapt

    ## compute the solution

    def calculateSolution(self,runName):
        """ Cacluate the PDEs numerical solution.

        Parameters
        ----------
        runName : str
            A name for the calculated solution.
        """
        logEvent("Setting initial conditions",level=0)
        for index,p,n,m,simOutput in zip(range(len(self.modelList)),self.pList,self.nList,self.modelList,self.simOutputList):
            if self.opts.hotStart:
                logEvent("Setting initial conditions from hot start file for "+p.name)
                tCount = int(self.ar[index].tree.getroot()[-1][-1][-1][0].attrib['Name'])
                self.ar[index].n_datasets = tCount + 1
                time = float(self.ar[index].tree.getroot()[-1][-1][-1][0].attrib['Value'])
                if len(self.ar[index].tree.getroot()[-1][-1]) > 1:
                    dt = time - float(self.ar[index].tree.getroot()[-1][-1][-2][0].attrib['Value'])
                else:
                    logEvent("Only one step in hot start file, setting dt to 1.0")
                    dt = 1.0
                logEvent("Last time step in hot start file was t = "+`time`)
                for lm,lu,lr in zip(m.levelModelList,m.uList,m.rList):
                    for cj in range(lm.coefficients.nc):
                        lm.u[cj].femSpace.readFunctionXdmf(self.ar[index],lm.u[cj],tCount)
                        lm.setFreeDOF(lu)
                        lm.timeIntegration.tLast = time
                        lm.timeIntegration.t = time
                        lm.timeIntegration.dt = dt
                self.tCount = tCount
            elif p.initialConditions is not None:
                logEvent("Setting initial conditions for "+p.name)
                m.setInitialConditions(p.initialConditions,self.tnList[0])
                #It's only safe to calculate the solution and solution
                #gradients because the models aren't attached yet
                for lm in m.levelModelList:
                    lm.calculateSolutionAtQuadrature()
            else:
                logEvent("No initial conditions provided for model "+p.name)
        if self.opts.hotStart:
            if time >= self.tnList[-1] - 1.0e-5:
                logEvent("Modifying time interval to be tnList[-1] + tnList since tnList hasn't been modified already")
                ndtout = len(self.tnList)
                dtout = (self.tnList[-1] - self.tnList[0])/float(ndtout-1)
                self.tnList = [time + i*dtout for i in range(ndtout)]
                logEvent("New tnList"+`self.tnList`)
            else:
                tnListNew=[time]
                for n,t in enumerate(self.tnList):
                    if time < t-1.0e-8:
                        tnListNew.append(t)
                self.tnList=tnListNew
                logEvent("Hotstarting, new tnList is"+`self.tnList`)
        else:
            self.tCount=0#time step counter
        logEvent("Attaching models and running spin-up step if requested")
        self.firstStep = True ##\todo get rid of firstStep flag in NumericalSolution if possible?
        spinup = []
        for index,m in self.modelSpinUp.iteritems():
            spinup.append((self.pList[index],self.nList[index],m,self.simOutputList[index]))
        for index,m in enumerate(self.modelList):
            logEvent("Attaching models to model "+p.name)
            m.attachModels(self.modelList)
            if index not in self.modelSpinUp:
                spinup.append((self.pList[index],self.nList[index],m,self.simOutputList[index]))
        for m in self.modelList:
            for lm,lu,lr in zip(m.levelModelList,
                                m.uList,
                                m.rList):
                #calculate the coefficients, any explicit-in-time
                #terms will be wrong
                lm.getResidual(lu,lr)
        for p,n,m,simOutput in spinup:
            logEvent("Attaching models to model "+p.name)
            m.attachModels(self.modelList)
            if m in self.modelSpinUp.values():
                logEvent("Spin-Up Estimating initial time derivative and initializing time history for model "+p.name)
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
                    lm.calculateAuxiliaryQuantitiesAfterStep()
                logEvent("Spin-Up Choosing initial time step for model "+p.name)
                m.stepController.initialize_dt_model(self.tnList[0],self.tnList[1])
                #mwf what if user wants spin-up to be over (t_0,t_1)?

                if m.stepController.stepExact and m.stepController.t_model_last != self.tnList[1]:
                    logEvent("Spin-up step exact called for model %s" % (m.name,),level=3)
                    m.stepController.stepExact_model(self.tnList[1])
                logEvent("Spin-Up Initializing time history for model step controller")
                m.stepController.initializeTimeHistory()
                m.stepController.setInitialGuess(m.uList,m.rList)
                solverFailed = m.solver.solveMultilevel(uList=m.uList,
                                                        rList=m.rList,
                                                        par_uList=m.par_uList,
                                                        par_rList=m.par_rList)
                Profiling.memory("solver.solveMultilevel")
                if solverFailed:
                    logEvent("Spin-Up Step Failed t=%12.5e, dt=%12.5e for model %s, CONTINUING ANYWAY!" %  (m.stepController.t_model,
                                                                                                     m.stepController.dt_model,
                                                                                                     m.name))

                else:
                    if n.restrictFineSolutionToAllMeshes:
                        logEvent("Using interpolant of fine mesh an all meshes")
                        self.restrictFromFineMesh(m)
                    self.postStep(m)
                    self.systemStepController.modelStepTaken(m,self.tnList[0])
                    logEvent("Spin-Up Step Taken, Model step t=%12.5e, dt=%12.5e for model %s" % (m.stepController.t_model,
                                                                                             m.stepController.dt_model,
                                                                                             m.name))
        for p,n,m,simOutput,index in zip(self.pList,self.nList,self.modelList,self.simOutputList,range(len(self.pList))):
            if not self.opts.hotStart:
                logEvent("Archiving initial conditions")
                self.archiveInitialSolution(m,index)
            else:
                self.ar[index].domain = self.ar[index].tree.find("Domain")
            self.initializeViewSolution(m)
            logEvent("Estimating initial time derivative and initializing time history for model "+p.name)
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
            logEvent("Choosing initial time step for model "+p.name)
            m.stepController.initialize_dt_model(self.tnList[0],self.tnList[1])
            #recalculate  with all terms ready
            for lm,lu,lr in zip(m.levelModelList,
                                m.uList,
                                m.rList):
                lm.getResidual(lu,lr)
            logEvent("Initializing time history for model step controller")
            m.stepController.initializeTimeHistory()
        self.systemStepController.initialize_dt_system(self.tnList[0],self.tnList[1]) #may reset other dt's
        for m in self.modelList:
            logEvent("Auxiliary variable calculations for model %s" % (m.name,))
            for av in self.auxiliaryVariables[m.name]:
                av.calculate_init()
        logEvent("Starting time stepping",level=0)
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

       # for p,n,m,simOutput,index in zip(self.pList,self.nList,self.modelList,self.simOutputList,range(len(self.pList))):
       #   for lm,lu,lr in zip(m.levelModelList,
       #                         m.uList,
       #                         m.rList):
       #     lm.getResidual(lu,lr)
       #     print "Initial Field \n %s" % lu
       #     print "Initial Residual \n %s" % lr
       #     print "Min / Max residual %s / %s" %(lr.min(),lr.max())

        self.nSequenceSteps = 0
        self.nSolveSteps=self.nList[0].adaptMesh_nSteps-1
        for (self.tn_last,self.tn) in zip(self.tnList[:-1],self.tnList[1:]):
            logEvent("==============================================================",level=0)
            logEvent("Solving over interval [%12.5e,%12.5e]" % (self.tn_last,self.tn),level=0)
            logEvent("==============================================================",level=0)
#            logEvent("NumericalAnalytics Time Step " + `self.tn`, level=0)
            if self.opts.save_dof:
                for m in self.modelList:
                    for lm in m.levelModelList:
                        for ci in range(lm.coefficients.nc):
                            lm.u[ci].dof_last[:] = lm.u[ci].dof
            if self.systemStepController.stepExact and self.systemStepController.t_system_last != self.tn:
                self.systemStepController.stepExact_system(self.tn)
            while self.systemStepController.t_system_last < self.tn:
                logEvent("System time step t=%12.5e, dt=%12.5e" % (self.systemStepController.t_system,
                                                              self.systemStepController.dt_system),level=3)

                while (not self.systemStepController.converged() and
                       not systemStepFailed):
                    logEvent("Split operator iteration %i" % (self.systemStepController.its,),level=3)
                    self.nSequenceSteps += 1
                    for (self.t_stepSequence,model) in self.systemStepController.stepSequence:
                        logEvent("NumericalAnalytics Model %s " % (model.name), level=0)
                        logEvent("Model: %s" % (model.name),level=1)
                        logEvent("NumericalAnalytics Time Step " + `self.t_stepSequence`, level=0)
                        logEvent("Fractional step %12.5e for model %s" % (self.t_stepSequence,model.name),level=3)
                        for m in model.levelModelList:
                            if m.movingDomain and m.tLast_mesh != self.systemStepController.t_system_last:
                                m.t_mesh = self.systemStepController.t_system_last
                                m.updateAfterMeshMotion()
                                m.tLast_mesh = m.t_mesh
                        self.preStep(model)
                        self.setWeakDirichletConditions(model)

                        stepFailed = False
                        if model.stepController.stepExact and model.stepController.t_model_last != self.t_stepSequence:
                            logEvent("Step exact called for model %s" % (model.name,),level=3)
                            model.stepController.stepExact_model(self.t_stepSequence)
                        while (model.stepController.t_model_last < self.t_stepSequence and
                               not stepFailed and
                               not self.systemStepController.exitModelStep[model]):
                            logEvent("Model step t=%12.5e, dt=%12.5e for model %s" % (model.stepController.t_model,
                                                                                 model.stepController.dt_model,
                                                                                 model.name),level=3)
                            for self.tSubstep in model.stepController.substeps:

                                logEvent("Model substep t=%12.5e for model %s" % (self.tSubstep,model.name),level=3)
                                #TODO: model.stepController.substeps doesn't seem to be updated after a solver failure unless model.stepController.stepExact is true
                                logEvent("Model substep t=%12.5e for model %s model.timeIntegration.t= %12.5e" % (self.tSubstep,model.name,model.levelModelList[-1].timeIntegration.t),level=3)

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
                                        logEvent("Using interpolant of fine mesh an all meshes")
                                        self.restrictFromFineMesh(model)
                                    model.stepController.updateSubstep()
                            #end model substeps
                            if solverFailed:
                                logEvent("Step failed due to solver failure")
                                stepFailed = not self.systemStepController.retryModelStep_solverFailure(model)
                            elif model.stepController.errorFailure():
                                logEvent("Step failed due to error failure")
                                stepFailed = not self.systemStepController.retryModelStep_errorFailure(model)
                            else:
                                #set up next step
                                self.systemStepController.modelStepTaken(model,self.t_stepSequence)                                
                                logEvent("Step Taken, t_stepSequence= %s Model step t=%12.5e, dt=%12.5e for model %s" % (self.t_stepSequence,
                                                                                                                         model.stepController.t_model,
                                                                                                                         model.stepController.dt_model,
                                                                                                                         model.name),level=3)
                        #end model step
                        if stepFailed:
                            logEvent("Sequence step failed")
                            if not self.systemStepController.ignoreSequenceStepFailure(model):
                                break
                            else:
                                logEvent("IGNORING STEP FAILURE")
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
                            logEvent("Retrying sequence")
                        else:
                            logEvent("Sequence failed")
                    else:
                        self.firstStep=False
                        systemStepFailed=False
                        logEvent("Step Taken, Model step t=%12.5e, dt=%12.5e for model %s" % (model.stepController.t_model,
                                                                                              model.stepController.dt_model,
                                                                                              model.name))
                        self.systemStepController.sequenceTaken()
                        for index,model in enumerate(self.modelList):
                            self.viewSolution(model,index)
                        if self.archiveFlag == ArchiveFlags.EVERY_MODEL_STEP:
                            self.tCount+=1
                            for index,model in enumerate(self.modelList):
                                self.archiveSolution(model,index,self.systemStepController.t_system)
                #end system split operator sequence
                if systemStepFailed:
                    logEvent("System Step Failed")
                    #go ahead and update as if the time step had succeeded
                    self.postStep(model)
                    self.systemStepController.modelStepTaken(model,self.t_stepSequence)
                    self.systemStepController.sequenceTaken()
                    self.systemStepController.updateTimeHistory()
                    #you're dead if retrySequence didn't work
                    logEvent("Step Failed, Model step t=%12.5e, dt=%12.5e for model %s" % (model.stepController.t_model,
                                                                                           model.stepController.dt_model,
                                                                                           model.name))
                    break
                else:
                    self.systemStepController.updateTimeHistory()
                    logEvent("Step Taken, System time step t=%12.5e, dt=%12.5e" % (self.systemStepController.t_system,
                                                                                   self.systemStepController.dt_system))
                    self.systemStepController.choose_dt_system()
                    logEvent("Potential System time step t=%12.5e, dt=%12.5e for next step" % (self.systemStepController.t_system,
                                                                                               self.systemStepController.dt_system))
                    if self.systemStepController.stepExact and self.systemStepController.t_system_last != self.tn:
                        self.systemStepController.stepExact_system(self.tn)
                for model in self.modelList:
                    for av in self.auxiliaryVariables[model.name]:
                        av.calculate()
                if self.archiveFlag == ArchiveFlags.EVERY_SEQUENCE_STEP:
                    self.tCount+=1
                    for index,model in enumerate(self.modelList):
                        self.archiveSolution(model,index,self.systemStepController.t_system_last)
                #can only handle PUMIDomain's for now
                self.nSolveSteps += 1
                if(self.PUMI_estimateError()):
                    self.PUMI_adaptMesh()
            #end system step iterations
            if self.archiveFlag == ArchiveFlags.EVERY_USER_STEP:
                self.tCount+=1
                for index,model in enumerate(self.modelList):
                    self.archiveSolution(model,index,self.systemStepController.t_system_last)
            if systemStepFailed:
                break
            #
            #h-adapt mesh, cekees modified from chitak
            #
            #assuming same for all physics and numerics  for now
        
            #can only handle PUMIDomain's for now
            self.nSolveSteps += 1
            if(self.PUMI_estimateError()):
              self.PUMI_adaptMesh()
        logEvent("Finished calculating solution",level=3)

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
            if (preCopy is not None and preCopy.has_key(('copy_uList')) and preCopy['copy_uList'] == True):
                for u_ci_lhs,u_ci_rhs in zip(levelModel.u.values(),self.modelList[preCopy['uList_model']].levelModelList[level].u.values()):
                    u_ci_lhs.dof[:] = u_ci_rhs.dof
                levelModel.setFreeDOF(model.uList[level])
            if preCopy is not None and preCopy.has_key(('clear_uList')) and preCopy['clear_uList'] == True:
                for u_ci_lhs in levelModel.u.values():
                    u_ci_lhs.dof[:] = 0.0
                levelModel.setFreeDOF(model.uList[level])
            if preCopy is not None and preCopy.has_key(('reset_uList')) and preCopy['reset_uList'] == True:
                levelModel.setFreeDOF(model.uList[level])
                levelModel.getResidual(model.uList[level],model.rList[level])

    ##intermodel transfer after a step
    def postStep(self,model):
        for level,levelModel in enumerate(model.levelModelList):
            postCopy = levelModel.coefficients.postStep(model.stepController.t_model,firstStep=self.firstStep)
            if postCopy is not None and postCopy.has_key(('copy_uList')) and postCopy['copy_uList'] == True:
                for u_ci_lhs,u_ci_rhs in zip(self.modelList[postCopy['uList_model']].levelModelList[level].u.values(),model.levelModelList[level].u.values()):
                    u_ci_lhs.dof[:] = u_ci_rhs.dof
                self.modelList[postCopy['uList_model']].levelModelList[level].setFreeDOF(self.modelList[postCopy['uList_model']].uList[level])

    def setWeakDirichletConditions(self,model):
        if model.weakDirichletConditions is not None:
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
        logEvent("Writing initial mesh for  model = "+model.name,level=3)
        logEvent("Writing initial conditions for  model = "+model.name,level=3)
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
        try:
            phi_s = {}
            phi_s[0] = model.levelModelList[-1].coefficients.phi_s
            model.levelModelList[-1].archiveFiniteElementResiduals(self.ar[index],
                                                                   self.tnList[0],
                                                                   self.tCount,
                                                                   phi_s,
                                                                   res_name_base='phi_s')
            logEvent("Writing initial phi_s at DOFs for = "+model.name+" at time t="+str(t),level=3)
        except:
            pass  

        #For aux quantity of interest (MQL)

        try:
            if model.levelModelList[-1].coefficients.outputQuantDOFs==True:
                try:
                    model.levelModelList[-1].setQuantDOFs()
                except:
                    pass
                quantDOFs = {}
                quantDOFs[0] = model.levelModelList[-1].quantDOFs
                model.levelModelList[-1].archiveFiniteElementResiduals(self.ar[index],
                                                                       self.tnList[0],
                                                                       self.tCount,
                                                                       quantDOFs,
                                                                       res_name_base='quantDOFs_for_'+model.name)
                logEvent("Writing initial quantity of interest at DOFs for = "+model.name+" at time t="+str(t),level=3)
        except:
            pass

        try:
            if model.levelModelList[-1].coefficients.outputQuantDOFs==True:
                try:
                    model.levelModelList[-1].setQuantDOFs()
                except:
                    pass
                quantDOFs2 = {}
                quantDOFs2[0] = model.levelModelList[-1].quantDOFs2
                model.levelModelList[-1].archiveFiniteElementResiduals(self.ar[index],
                                                                       self.tnList[0],
                                                                       self.tCount,
                                                                       quantDOFs2,
                                                                       res_name_base='quantDOFs2_for_'+model.name)
                logEvent("Writing initial quantity of interest at DOFs for = "+model.name+" at time t="+str(t),level=3)
        except:
            pass           

        #Write bathymetry for Shallow water equations (MQL)
        try:
            bathymetry = {}
            bathymetry[0] = model.levelModelList[-1].coefficients.b.dof
            model.levelModelList[-1].archiveFiniteElementResiduals(self.ar[index],
                                                                   self.tnList[0],
                                                                   self.tCount,
                                                                   bathymetry,
                                                                   res_name_base='bathymetry')
            logEvent("Writing bathymetry for = "+model.name,level=3)
        except:
            pass
        #write eta=h+bathymetry for SWEs (MQL)
        try:
            eta = {}
            eta[0] = model.levelModelList[-1].coefficients.b.dof+model.levelModelList[-1].u[0].dof
            model.levelModelList[-1].archiveFiniteElementResiduals(self.ar[index],
                                                                   self.tnList[0],
                                                                   self.tCount,
                                                                   eta,
                                                                   res_name_base='eta')
            logEvent("Writing bathymetry for = "+model.name,level=3)
        except:
            pass

        #for nonlinear POD
        if self.archive_pod_residuals[index] == True:
            res_space = {}; res_mass = {}
            for ci in range(model.levelModelList[-1].coefficients.nc):
                res_space[ci] = numpy.zeros(model.levelModelList[-1].u[ci].dof.shape,'d')
                model.levelModelList[-1].getSpatialResidual(model.levelModelList[-1].u[ci].dof,res_space[ci])
                res_mass[ci] = numpy.zeros(model.levelModelList[-1].u[ci].dof.shape,'d')
                model.levelModelList[-1].getMassResidual(model.levelModelList[-1].u[ci].dof,res_mass[ci])
            model.levelModelList[-1].archiveFiniteElementResiduals(self.ar[index],self.tnList[0],self.tCount,res_space,res_name_base='spatial_residual')
            model.levelModelList[-1].archiveFiniteElementResiduals(self.ar[index],self.tnList[0],self.tCount,res_mass,res_name_base='mass_residual')

        if not self.opts.cacheArchive:
            if not self.so.useOneArchive:
                self.ar[index].sync()
            else:
                if index == len(self.ar) - 1:
                    self.ar[index].sync()
    ##save model's solution values to archive
    def archiveSolution(self,model,index,t=None):
        if self.archiveFlag == ArchiveFlags.UNDEFINED:
            return
        if t is None:
            t = self.systemStepController.t_system

        logEvent("Writing mesh header for  model = "+model.name+" at time t="+str(t),level=3)
        logEvent("Writing solution for  model = "+model.name,level=3)
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

        try:
            phi_s = {}
            phi_s[0] = model.levelModelList[-1].coefficients.phi_s
            model.levelModelList[-1].archiveFiniteElementResiduals(self.ar[index],
                                                                   self.tnList[0],
                                                                   self.tCount,
                                                                   phi_s,
                                                                   res_name_base='phi_s')
            logEvent("Writing phi_s at DOFs for = "+model.name+" at time t="+str(t),level=3)
        except:
            pass
        
        try:
            if model.levelModelList[-1].coefficients.outputQuantDOFs==True:
                quantDOFs = {}
                quantDOFs[0] = model.levelModelList[-1].quantDOFs
                model.levelModelList[-1].archiveFiniteElementResiduals(self.ar[index],
                                                                       self.tnList[0],
                                                                       self.tCount,
                                                                       quantDOFs,
                                                                       res_name_base='quantDOFs_for_'+model.name)        
                logEvent("Writing quantity of interest at DOFs for = "+model.name+" at time t="+str(t),level=3)
        except:
            pass

        try:
            if model.levelModelList[-1].coefficients.outputQuantDOFs==True:
                quantDOFs2 = {}
                quantDOFs2[0] = model.levelModelList[-1].quantDOFs2
                model.levelModelList[-1].archiveFiniteElementResiduals(self.ar[index],
                                                                       self.tnList[0],
                                                                       self.tCount,
                                                                       quantDOFs2,
                                                                       res_name_base='quantDOFs2_for_'+model.name)        
                logEvent("Writing quantity of interest at DOFs for = "+model.name+" at time t="+str(t),level=3)
        except:
            pass        

        #Write bathymetry for Shallow water equations (MQL)
        try:
            bathymetry = {}
            bathymetry[0] = model.levelModelList[-1].coefficients.b.dof
            model.levelModelList[-1].archiveFiniteElementResiduals(self.ar[index],
                                                                   self.tnList[0],
                                                                   self.tCount,
                                                                   bathymetry,
                                                                   res_name_base='bathymetry')
            logEvent("Writing bathymetry for = "+model.name,level=3)
        except:
            pass
        #write eta=h+bathymetry for SWEs (MQL)
        try:
            eta = {}
            eta[0] = model.levelModelList[-1].coefficients.b.dof+model.levelModelList[-1].u[0].dof
            model.levelModelList[-1].archiveFiniteElementResiduals(self.ar[index],
                                                                   self.tnList[0],
                                                                   self.tCount,
                                                                   eta,
                                                                   res_name_base='eta')
            logEvent("Writing bathymetry for = "+model.name,level=3)
        except:
            pass

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
            if not self.so.useOneArchive:
                self.ar[index].sync()
            else:
                if index == len(self.ar) - 1:
                    self.ar[index].sync()

    ## clean up archive
    def closeArchive(self,model,index):
        if self.archiveFlag is None:
            return
        if self.so.useOneArchive:
            if index==0:
                logEvent("Closing solution archive for "+self.so.name)
                self.ar[index].close()
        else:
            logEvent("Closing solution archive for "+model.name)
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
