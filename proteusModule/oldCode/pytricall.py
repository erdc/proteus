#!/usr/bin/env python
#import standard Python modules
import sys,os
#module for triangle
import triangulate
import numpy as Numeric



def tricall0():
    import triIfaceUtils
    import triIfaceFileUtils
    #define some flags for how the script will run
    verbose = 2   #level of output
    tinfo = """triangulation:
number of points        = %d
attributes per point    = %d
number of triangles     = %d
nodes per triangle      = %d
attributes per triangle = %d
number of segments      = %d
number of holes         = %d
number of regions       = %d
number of edges         = %d
"""

    trin = triangulate.new()
    pointsin = Numeric.array([[0.0, 0.0],
                              [1.0,0.0],
                              [1.0,10.0],
                              [0.0,10.0]])
    npoints = pointsin.shape[0]
    pattsin = Numeric.array([0.0, 1.0, 11.0, 10.0])
    pmarkin = Numeric.zeros((npoints,),Numeric.Int)
    pmarkin[1] = 2
    pregion =  Numeric.array([[0.5,5.0,7.0,0.1]])
    if verbose > 3:
        print 'pregion shape=',pregion.shape
    #end verbose
    triangulate.setPointsAndMarkers(trin,pointsin,pmarkin)
    triangulate.setPointAttributes(trin,pattsin)
    triangulate.setRegions(trin,pregion)
    
    tiInfo = triangulate.getInfo(trin)

    print 'trin info says'
    print tinfo % tiInfo

    #create a simple flag list
    flags = "pczAevn"
    
    trimid = triangulate.new()
    vorout = triangulate.new()
    
    triangulate.applyTriangulate(flags,trin,trimid,vorout)

    tmInfo = triangulate.getInfo(trimid)

    print 'trimid info says'
    print tinfo % tmInfo

    if verbose > 1:
        points0 = triangulate.getPoints(trin)
        print 'in points=\n',points0
        elems0 = triangulate.getTriangles(trin)
        print 'in elems=\n',elems0
        elems1 = triangulate.getTriangles(trimid)
        print 'out elems=\n',elems1
    #end verbose
    if verbose > 1:
        print 'calling printReport for input'
        triangulate.printReport(trin)
        print 'calling printReport for output'
        triangulate.printReport(trimid)
    #end verbose
    triIfaceUtils.writeOutTriangulation(trimid,"trimesh")

    #now refine intermediate triangle
    nelmsmid = tmInfo[2]
    midArea = Numeric.zeros((nelmsmid,),Numeric.Float)
    midArea[0] = 3.0
    midArea[1] = 1.0

    triangulate.setTriangleAreas(trimid,midArea)
    
    triout = triangulate.new()
    flags  = "prazBP"
    triangulate.applyTriangulateNoVoronoi(flags,trimid,triout)
    
    triIfaceUtils.writeOutTriangulation(triout,"trimesh.1")


    #now view with showme
    showme = "../bin/showme"
    showmecmd = showme+" trimesh"

    failure = 0
    failure = triIfaceFileUtils.checkFileExists(showme)
    failure = os.system(showmecmd)
    
    #manually delete things for debugging
    print 'deleting trin'
    del trin
    print 'deleting trimid'
    del trimid
    print 'deleting vor'
    del vorout
    print 'deleting triout'
    del triout

    return failure
#end tricall0


def tricall2():
    """
    test reading Triangle file utilities
    """
    import triIfaceUtils
    import triIfaceFileUtils

    reader = triIfaceUtils.TriangleInputFileReader()
    filebase = 'trimesh'
    
    nodeDataInfo,nodeData = reader.readNodes(filebase)
    nodes,nodesA,nodesM = (nodeData['nodes'],
                           nodeData['nodeAttributes'],
                           nodeData['nodeMarkers'])
    
    print 'Nodes: nodeInfo= ',nodeDataInfo
    print """Nodes: nodes= \n%s\n nodeAttributes= \n%s\n nodeMarkers= \n%s\n """ % (nodes,
                                                                                    nodesA,
                                                                                    nodesM)
    
    triDataInfo,triData = reader.readTriangles(filebase)
    triangles,trianglesA = triData['triangles'],triData['triangleAttributes']
    
    print 'Triangles: triInfo= ',triDataInfo
    print """Triangles: elems= \n%s\n triAttributes= \n%s\n""" % (triangles,
                                                                  trianglesA)

    #write out assuming just read in node and ele files
    filebaseout = filebase+'-out-0'
    tri0 = triangulate.new()

    
    if nodesM == None:
        triangulate.setPoints(tri0,nodes)
    else:
        triangulate.setPointsAndMarkers(tri0,nodes,nodesM)
    if not nodesA == None: 
        triangulate.setPointAttributes(tri0,nodesA)
    #end if
    
    triangulate.setTriangles(tri0,triangles)
    if not trianglesA == None: 
        triangulate.setTriangleAttributes(tri0,trianglesA)
    #end if
    triIfaceUtils.writeOutTriangulation(tri0,filebaseout,nbase=0,verbose=0)

    #now view with showme
    showme = "../bin/showme"
    showmecmd = """showme %s """ % filebaseout

    failure = 0
    failure = triIfaceFileUtils.checkFileExists(showme)
    failure = os.system(showmecmd)
    
    ##now read from .poly file
    polyDataInfo,polyData = reader.readPoly(filebase)

    
    for type in ['node','segment','hole','region']:
        print 'Poly: ',type,'Info= \n',polyDataInfo[type],'\n'
    #
    nodes2,nodesA2,nodesM2 = (polyData['node']['nodes'],
                              polyData['node']['nodeAttributes'],
                              polyData['node']['nodeMarkers'])
    segments2,segmentsM2   = (polyData['segment']['segments'],
                              polyData['segment']['segmentMarkers'])
    holes2                 = polyData['hole']['holes']
    regions2               = polyData['region']['regions']
    
    print """Poly file read:
nodes   = \n%s\n nodeAttributes= \n%s\n nodeMarkers= \n%s\n 
segments= \n%s\n segmentMarkers= \n%s\n 
holes   = \n%s\n
regions = \n%s\n
    """ % (nodes2,nodesA2,nodesM2,segments2,segmentsM2,holes2,regions2)

    ### now create a triangulation and write out the data arrays?
    trin1 = triangulate.new()
    
    filebaseout1 = filebase+"-out-1"

    if nodesM2 == None:
        triangulate.setPoints(trin1,nodes2,nodesM2)
    else:
        triangulate.setPointsAndMarkers(trin1,nodes2,nodesM2)
    if not nodesA2 == None:
        triangulate.setPointAttributes(trin1,nodesA2)
    #end if

    if segmentsM2 == None:
        triangulate.setSegments(trin1,segments2)
    else:
        triangulate.setSegmentsAndMarkers(trin1,segments2,segmentsM2)
    if not holes2 == None:
        triangulate.setHoles(trin1,holes2)
    #end if
    if not regions2 == None:
        #print 'setting trin1 regions=\n',regions2
        triangulate.setRegions(trin1,regions2)
    #end if
    trout1 = triangulate.new()
    flags = "zpne" #needs to be base 0 if Input was that way (and nbase=0)
    triangulate.applyTriangulateNoVoronoi(flags,trin1,trout1)

    triIfaceUtils.writeOutTriangulation(trout1,filebaseout1,nbase=0,verbose=0)

    #now view with showme
    showmecmd = """showme %s """ % filebaseout1

    failure = 0
    failure = triIfaceFileUtils.checkFileExists(showme)
    failure = os.system(showmecmd)

#end tricall2

def tricall3():
    import TriangleIface
    verbose = 1
    filebase = "trimesh"
    mesh = TriangleIface.TriangleBaseMesh(verbose=verbose)

    #read in from .node and .ele files
    flagsAdd= "" # "r"
    
    mesh.readFromNodeAndEleFiles(filebase,flagsAdd,verbose=verbose)

    print 'viewing mesh generated from .node and .ele files'
    mesh.viewShowme()
    
    mesh1 = TriangleIface.TriangleBaseMesh(verbose=verbose)

    #read in from .node and .ele files
    flagsAdd= "p" # "p"
    
    mesh1.readFromPolyFile(filebase,flagsAdd,verbose=verbose)

    print 'viewing mesh1 generated from .poly file'
    mesh1.viewShowme()
#end triIfaceCall3
def exPyadhLaplace1(filebase="trimesh",baseFlags="zen",
                    flagsAdd="",viewMesh=1,verbose=0):
    import QuadTools
    import FemTools
    import TriangleIface
    import PoissonTestProblems
    import ScalarTransport
    import TimeIntegrationTools
    import numpy as numpy.oldnumeric.linear_algebra as LinearAlgebraTools
    import LinearSolvers
    import NonlinearSolvers
    
    nbase = 0
    if baseFlags.find('z') == -1:
        nbase=1
    mesh = TriangleIface.TriangleBaseMesh(baseFlags=baseFlags,
                                          nbase=nbase,
                                          verbose=verbose)

    if flagsAdd.find('p') >= 0:
        if verbose > 0:
            print """flagsAdd= %s trying to read %s """ % (flagsAdd,
                                                           filebase+'.poly')
        #end verbose
        mesh.readFromPolyFile(filebase,flagsAdd,verbose=verbose)

    elif flagsAdd.find('r') >= 0:
        if verbose > 0:
            print """flagsAdd= %s trying to read %s and %s """ % (flagsAdd,
                                                                  filebase+'.node',
                                                                  filebase+'.ele')
        #end verbose
        mesh.readFromNodeAndEleFiles(filebase,flagsAdd,verbose=verbose)
        
    else:
        if verbose > 0:
            print """flagsAdd= %s trying to read %s """ % (flagsAdd,
                                                           filebase+'.node')
        #end verbose
        
        mesh.readFromNodeFile(filebase,flagsAdd,verbose=verbose)
    #end if on flags
    if viewMesh > 0:
        print 'viewing mesh generated from .node and .ele files'
        mesh.viewShowme()
    #end if
    pyadhMesh = mesh.convertToPyadhMesh(verbose)
    pyadhMesh.writeEdgesGnuplot2("gnuMesh") #uses array interface
    if viewMesh > 0:
        pyadhMesh.viewMeshGnuplotPipe("gnuMesh")
    matmesh = "matlabMesh"
    pyadhMesh.buildMatlabMeshDataStructures(matmesh)
    pyadhMesh.writeMeshEnsight("ex1soln")
    #####
    #solve a simple poisson equation with Dirichlet bc's on left and right?
    nd = 2
    f = lambda p : 0.0
    def dirBCs(p):
        if p[0] == -1.0:
            return lambda p,t: 10.0
        elif p[0] >= 40.8893:
            return lambda p,t: 0.0
    #end dirBCs
    A  = Numeric.zeros((2,2),Numeric.Float)
    A[0,0] = 1.0; A[1,1]= 1.0
    pc = PoissonTestProblems.PoissonCoefficients(A,f,nd)
    dirProb = dirBCs
    useCG = True
    order = 1
    quadratureOrder=3
    bquadratureOrder=3
    quadrature = {}
    #end quadrature order
    elemQuadrature = QuadTools.SimplexGaussQuadrature(nd)
    elemQuadrature.setOrder(quadratureOrder)
    for integral in ScalarTransport.OneLevelScalarTransport.integralKeys:
        quadrature[integral] = elemQuadrature
    #end for
    elementBoundaryQuadrature = {}
    boundaryQuadrature = QuadTools.SimplexGaussQuadrature(nd-1)
    boundaryQuadrature.setOrder(bquadratureOrder)
    for bintegral in ScalarTransport.OneLevelScalarTransport.elementBoundaryIntegralKeys:
        elementBoundaryQuadrature[bintegral] = boundaryQuadrature
    ### setup finite element spaces and stuff
    if verbose > 0:
        print 'definining finite element spaces'
    #end verbose
    #try P^2
    if useCG:
        if order == 1:
            FemSpace = FemTools.C0_AffineLinearOnSimplexWithNodalBasis(pyadhMesh,nd)
        else:
            FemSpace = FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis(pyadhMesh,nd)
    else:
        if order == 1:
            FemSpace = FemTools.DG_AffineLinearOnSimplexWithNodalBasis(pyadhMesh,nd)
        else:
            FemSpace = FemTools.DG_AffineQuadraticOnSimplexWithNodalBasis(pyadhMesh,nd)

    #end if
    u   = FemTools.FiniteElementFunction(FemSpace)
    phi = FemTools.FiniteElementFunction(FemSpace)
    #phi = u
    if verbose > 0:
        print 'definining boundary conditions'
    #end verbose
    dirichletBCs=FemTools.DOFBoundaryConditions(
        FemSpace,dirProb)
    fluxBndyCond='noFlow'
    ### setup numerical approximation configuration
    if verbose > 0:
        print 'definining linear algebra'
    #end verbose
    MatType = LinearAlgebraTools.SparseMat
    matType = 'csr'
    linearSolverType=  'SparseLU'
    ### time integration
    #steady state
    if verbose > 0:
        print 'definining time integration'
    #end verbose
    TimeIntegrationClass = TimeIntegrationTools.TimeIntegration
    ### flux approximations
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
    ###create a single level solver
    if verbose > 0:
        print 'creating single level system'
    #end verbose
    system = ScalarTransport.OneLevelScalarTransport(u,
                                                     phi,
                                                     FemSpace,
                                                     dirichletBCs,
                                                     pc,
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
    
    ### create nonlinear system for solving problem
    if verbose > 0:
        print 'creating nonlinear system'
    #end verbose
    
    jacobian = MatType(system.dim,system.dim,min(7,system.dim))
    #mwf add following Chris example
    jacobian = system.initializeJacobian(jacobian)
    ### create linear system for solving problem
    if verbose > 5:
        print 'creating linear system'
        print 'initial jacobian mat = ',jacobian
    #end verbose

    linearSolver = LinearSolvers.SparseLU(jacobian)

    nlSolver = NonlinearSolvers.Newton(linearSolver,system,jacobian,maxIts=10)
    #need intial residual calculated
    system.getResidual(y,r)
    if verbose > 5:
        print 'system y0= ',y
        system.getJacobian(jacobian)
        print 'system residual= ',r
        print 'system jacobian = ',jacobian
    ### solving problem
    if verbose > 0:
        print 'trying to solve nonlinear system'
    #end verbose
    FemSpace.writeFunctionEnsight(u,"ex1soln",append=False)
    nlSolver.solve(y,r)

        
    if verbose > 0:
        print 'nonlinear solver done'
    #end if

    #print out solution
    FemSpace.writeFunctionMatlab(u,"ex1soln")
    FemSpace.writeFunctionEnsight(u,"ex1soln",append=True)
    FemSpace.endTimeSeriesEnsight([0,1],"ex1soln","ex1soln")
#end exAdh
    
if __name__ == '__main__':
    #make sure python has pyadh in path
    import os,sys
    from optparse import OptionParser
    parser = OptionParser()

    #options controlling simulation behavior
    parser.add_option('-P','--pyadhDir',
                      default='/Users/mfarthin/Public/code/chrisAdhPyUtil/cek-dev/',
                      help="""where to find pyadh library""")

    parser.add_option('-v','--verbose',
                      default=0,
                      help="""level of verbosity in simulation [0]""")

    #get options
    options, args = parser.parse_args(sys.argv[1:]) #get command line args

    #
    verbose  = int(options.verbose)
    pyadhDir = str(options.pyadhDir)
    sys.path.insert(0,pyadhDir)

    #tricall0()
    #tricall2()
    #tricall3()
    flagsAdd="pqa" #"pq"
    
    exPyadhLaplace1(filebase=os.getenv('PYADH_PACKAGES',
                                           os.getenv('HOME')+
                                           '/src/pyadh-packages')+
                        "/triangle/examples/la",
                        baseFlags="en",
                        flagsAdd=flagsAdd,viewMesh=1,verbose=0)
    
