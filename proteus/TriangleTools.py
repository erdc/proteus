#!/usr/bin/env python
"""Collect functionality to make a lightweight interface for
Shewchuk's Triangle package. Follow the ellipt2d example in some
respects

.. inheritance-diagram:: proteus.TriangleTools
   :parts: 1
"""
import os
import numpy
import triangleWrappers
import TriangleUtils
import TriangleFileUtils
import MeshTools

def genMeshWithTetgen(polyfile,
                      baseFlags="Yp",
                      nbase = 1):
    """
    Generate a tetrahedral mesh from a polyfile.

    Arguments
    ---------
    polyfile : str
        Filename with appropriate data for tengen.
    baseFlags : str
        Standard Tetgen options for generation
    nbase : int

    Returns
    -------
    mesh : :class:`proteus.MeshTools.TetrahedralMesh`
        Simplex mesh
    """
    from subprocess import check_call
    tetcmd = "tetgen - %s %s.poly" % (baseFlags, polyfile)
    check_call(tetcmd,shell=True)
    elefile = "%s.1.ele" % polyfile
    nodefile = "%s.1.node" % polyfile
    facefile = "%s.1.face" % polyfile
    edgefile = "%s.1.edge" % polyfile
    assert os.path.exists(elefile), "no 1.ele"
    tmp = "%s.ele" % polyfile
    os.rename(elefile,tmp)
    assert os.path.exists(tmp), "no .ele"
    assert os.path.exists(nodefile), "no 1.node"
    tmp = "%s.node" % polyfile
    os.rename(nodefile,tmp)
    assert os.path.exists(tmp), "no .node"
    if os.path.exists(facefile):
        tmp = "%s.face" % polyfile
        os.rename(facefile,tmp)
        assert os.path.exists(tmp), "no .face"
    if os.path.exists(edgefile):
        tmp = "%s.edge" % polyfile
        os.rename(edgefile,tmp)
        assert os.path.exists(tmp), "no .edge"
    mesh=MeshTools.TetrahedralMesh()
    mesh.generateFromTetgenFiles(polyfile,
                                 base=nbase)
    return mesh



class TriangleBaseMesh:
    """ This is basically a wrapper for the triangulation interface
    that should be able to create a triangle mesh in different ways

       from .ele and .node files
       from a .poly file
       from a proteus mesh

    It should also be able to generate a proteus mesh from the
    triangle representation and allow the user to access the basic
    data arrays in triangle

    """

    def __init__(self,baseFlags="zen",nbase=0,verbose=0):
        """
        initialize the triangulation object,
        keep track of what numbering scheme it uses,
        create a base set of flags for triangulate (e.g., z if using base 0)

        does not create a Voronoi diagram by default
        """
        self.trirep = []
        self.trirep.append(triangleWrappers.new())
        self.trirepDefaultInit = []
        self.trirepDefaultInit.append(True)
        self.vorrep = []
        self.vorrep.append(triangleWrappers.new()) #don't create a Voronoi diagram by default
        self.makeVoronoi = False
        self.vorrepDefaultInit = []
        self.vorrepDefaultInit.append(True)
        assert(nbase in [0,1])
        self.nbase  = nbase
        self.baseFlags = baseFlags
        if self.nbase == 0 and self.baseFlags.find('z') == -1:
            print """WARNING TriangleMesh base numbering= %d
            reqires z in baseFlags = %s, adding """ % (self.nbase,self.baseFlags)
            self.baseFlags += "z"
        #end if
        if self.baseFlags.find('v') >= 0:
            self.makeVoronoi = True
        if verbose > 0:
            print """TriangleBaseMesh nbase=%d baseFlags= %s """ % (self.nbase,
                                                                    self.baseFlags)
        #end if
        #keep track of last set of flags used to generate rep
        self.lastFlags = None
        #to make printing out info easier
        self.infoFormat = """triangulation:
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

        #end init
    def resetDefaultTrirep(self,index=0,verbose=0):
        """
        reset the Triangle mesh representation if it has been set to something nontrivial

        """
        assert(index in range(0,len(self.trirep)))
        #delete the existing representation if it has been initialized
        if not self.trirepDefaultInit[index]:
            if verbose > -1:
                print 'TriangleIface deleting current trirep[',index,']'
            #end verbose
            del self.trirep[index]
            #get a new represenation
            if index == 0:
                self.trirep.append(triangleWrappers.new())
            else:
                self.trirep[index] = triangleWrappers.new()
            self.trirepDefaultInit[index] = True
        #end trirep initialized to something non trivial

    def resetDefaultVorrep(self,index=0,verbose=0):
        """
        reset the Triangle mesh Voronoi representation if it has
        been set to something nontrivial
        """
        assert(index in range(0,len(self.vorrep)))
        #delete the existing representation if it has been initialized
        if not self.vorrepDefaultInit[index]:
            if verbose > -1:
                print 'TriangleIface deleting current vorrep[',index,']'
            #end verbose
            del self.vorrep[index]
            #get a new represenation
            self.vorrep[index] = triangleWrappers.new()
            self.vorrepDefaultInit[index] = True
        #end trirep initialized to something non trivial

    def convertToProteusMesh(self,verbose=0):
        """
        Generate a representation in the format expected by
        proteus.

        Need to make sure the triangle mesh has generated
           nodes
           elements (triangles)
           edges
           neighbors

        First set the _global arrays for
            nodes
            elements
        """
        import MeshTools
        triInfo = triangleWrappers.getInfo(self.trirep[0])
        if verbose > 1:
            print 'generating proteus Mesh:'
            print self.infoFormat % triInfo
        #end if

        #create basic proteus mesh
        meshout = MeshTools.TriangularMesh()
        meshout.generateFromTriangleMesh(self.trirep[0],self.nbase)
        return meshout

#         #get basic information to make sure I can proceed
#         #with current mesh representation
#         nNodes_global = triInfo[0]; nElems_global = triInfo[2];
#         nNodes_elem   = triInfo[3]; nEdges_global = triInfo[8];
#         spaceDim      = 2
#         assert(nNodes_global > 0 and nElems_global > 0
#                and nNodes_elem >= 3 and nEdges_global > 0)

#         #subtract off base since proteus.Mesh wants base 0 more or less
#         nbase = self.nbase
#         #should also check? base zero,
#         #get the minimum array information
#         ##nodes
#         nodeArray = triangleWrappers.getPoints(self.trirep[0])
#         meshout.nNodes_global     = nNodes_global
#         meshout.nodeArray = numpy.zeros((meshout.nNodes_global,3),
#                                           'd')
#         for nN in range(nNodes_global):
#             meshout.nodeArray[nN,:spaceDim] = nodeArray[nN,:spaceDim]
#             #cek added because refinement still depends on list interface
#             #mwf commented  out 01/29/08
#             #meshout.nodeList.append(MeshTools.Node(nN,meshout.nodeArray[nN,0],meshout.nodeArray[nN,1]))
#         ##elements
#         #cek added nodeDict init b/c of refinement
#         #mwf commented out 01/29/08
#         #meshout.nodeDict = dict([(n,n) for n in meshout.nodeList])
#         elemArray  = triangleWrappers.getTriangles(self.trirep[0])
#         #ignore higher order nodes for proteus mesh
#         meshout.nNodes_element    = spaceDim+1
#         meshout.nElements_global  = nElems_global
#         #just copy over first 3 nodes
#         meshout.elementNodesArray = numpy.zeros((nElems_global,spaceDim+1),
#                                                   'i')
#         for eN in range(nElems_global):
#             meshout.elementNodesArray[eN,:] = elemArray[eN,:spaceDim+1]-nbase
#             #cek added b/c of refinement
#             #mwf commented out 01/29/08
#             #meshout.newTriangle([meshout.nodeList[meshout.elementNodesArray[eN,0]],
#             #                     meshout.nodeList[meshout.elementNodesArray[eN,1]],
#             #                     meshout.nodeList[meshout.elementNodesArray[eN,2]]])
#         #end eN
#         #cek added finalize call to build lists and arrays
#         #mwf commented out 01/29/08
#         #meshout.finalize()
#         ##proteusMesh keeps elements per node
# #         nodeElementsDict={}
# #         for eN in range(nElems_global):
# #             for nN_element in range(spaceDim+1):
# #                 nN = meshout.elementNodesArray[eN,nN_element]
# #                 if  nodeElementsDict.has_key(nN):
# #                     nodeElementsDict[nN].append(eN)
# #                 else:
# #                     nodeElementsDict[nN] = [eN]
# #                 #end if
# #             #end for nN_element
# #         #end eN
# #         meshout.max_nElements_node = max(len(nodeElementsDict[nN]) for nN in range(meshout.nNodes_global))

# #         meshout.nElements_node    = numpy.zeros((meshout.nNodes_global,),'i')
# #         meshout.nodeElementsArray = numpy.zeros((meshout.nNodes_global,
# #                                                    meshout.max_nElements_node),
# #                                                   'i')
# #         for nN,elementList in nodeElementsDict.iteritems():
# #             meshout.nElements_node[nN] = len(elementList)
# #             for eN_element,eN in enumerate(elementList):
# #                 meshout.nodeElementsArray[nN,eN_element]=eN
# #             #end eN_element
# #         #end nN

# #         ##now build the element <--> element boundary arrays
# #         meshout.nElementBoundaries_element = spaceDim+1
# #         #make sure Triangle keeps all edges and not just boundary ones
# #         meshout.nElementBoundaries_global  = nEdges_global

# #         #maps element, local edge number to global edge number
# #         meshout.elementBoundariesArray = numpy.zeros((nElems_global,spaceDim+1),
# #                                                        'i')

# #         #maps edge to global element on left and right (out of domain is 0)
# #         meshout.elementBoundaryElementsArray=numpy.ones(
# #             (meshout.nElementBoundaries_global,2),'i')
# #         meshout.elementBoundaryElementsArray*=-1
# #         #maps global edge to its local number on the neighboring elements
# #         meshout.elementBoundaryLocalElementBoundariesArray=numpy.zeros(
# #             (meshout.nElementBoundaries_global,2),'i')

# #         #several options for edge to "left" and "right" element neighbor,
# #         #could number each neighbor according to which one
# #         #has a given edge first in the current numbering
# #         #could try to make element 0 be the one that has same orientation of edge

# #         elementBoundaryElementsCardArray = numpy.zeros((meshout.nElementBoundaries_global,),
# #                                                          'i') #CardArray in Mesh
# #         #I have to generate the element-->global edge number mapping manually
# #         edgeArray = triangleWrappers.getEdges(self.trirep[0])
# #         edgeDict  = {}
# #         for edgeN in range(nEdges_global):
# #             n0 = edgeArray[edgeN,0]-nbase; n1 = edgeArray[edgeN,1]-nbase
# #             edgeDict[(n0,n1)] = edgeN #global edge number
# #         #end for

# #         for eN in range(nElems_global):
# #             locNodes = elemArray[eN,:spaceDim+1]-nbase #global node numbers
# #             edgeOrientedSame = [False,False,False]
# #             #edge I is across from node I
# #             #0
# #             e0 = (locNodes[1],locNodes[2]); e0rev = (locNodes[2],locNodes[1])
# #             e0_global = None
# #             if edgeDict.has_key(e0):
# #                 e0_global = edgeDict[e0] #same orientation
# #                 edgeOrientedSame[0] = True
# #             elif edgeDict.has_key(e0rev):     #opposite orientation
# #                 e0_global = edgeDict[e0rev]   #could make negative to denote orientation
# #             #end if
# #             assert(not e0_global == None)
# #             #1
# #             e1 = (locNodes[2],locNodes[0]); e1rev = (locNodes[0],locNodes[2])
# #             e1_global = None
# #             if edgeDict.has_key(e1):
# #                 e1_global = edgeDict[e1] #same orientation
# #                 edgeOrientedSame[1] = True
# #             elif edgeDict.has_key(e1rev):     #opposite orientation
# #                 e1_global = edgeDict[e1rev]   #could make negative to denote orientation
# #             #end if
# #             assert(not e1_global == None)
# #             #2
# #             e2 = (locNodes[0],locNodes[1]); e2rev = (locNodes[1],locNodes[0])
# #             e2_global = None
# #             if edgeDict.has_key(e2):
# #                 e2_global = edgeDict[e2] #same orientation
# #                 edgeOrientedSame[2] = True
# #             elif edgeDict.has_key(e2rev):     #opposite orientation
# #                 e2_global = edgeDict[e2rev]   #could make negative to denote orientation
# #             #end if
# #             assert(not e2_global == None)
# #             eI_global = numpy.array([e0_global,e1_global,e2_global],
# #                                       'i')
# #             meshout.elementBoundariesArray[eN,:] = eI_global
# #             for eI in eI_global:
# #                 #edge --> element mappings
# #                 elementBoundaryElementsCardArray[eI] += 1
# #             #end eI


# #             for eiloc in range(meshout.nElementBoundaries_element):
# #                 eI     = eI_global[eiloc]
# #                 #first visited labelling
# #                 elneig = elementBoundaryElementsCardArray[eI]-1
# #                 #elneig = 0 #same orientation
# #                 #if not edgeOrientedSame[eiloc]
# #                 #elneig = 1 #opposite
# #                 ##endif
# #                 meshout.elementBoundaryElementsArray[eI,elneig]=eN
# #                 meshout.elementBoundaryLocalElementBoundariesArray[eI,elneig]=eiloc
# #             #end eiloc
# #         #end eN

# #         #perform some sanity checks
# #         for edgeN in range(nEdges_global):
# #             assert(elementBoundaryElementsCardArray[edgeN] in [1,2])
# #             eN0 = meshout.elementBoundaryElementsArray[edgeN,0]
# #             eN1 = meshout.elementBoundaryElementsArray[edgeN,1]
# #             assert(0 <= eN0 and eN0 < meshout.nElements_global)
# #             if elementBoundaryElementsCardArray[edgeN] == 1:
# #                 assert(eN1 == -1)
# #             else:
# #                 assert(0 <= eN1 and eN1 < meshout.nElements_global)
# #             #end if
# #         #end sanity check on edges
# #         #now take care of boundary edges
# #         #interior elements counted twice
# #         sumCard = numpy.sum(elementBoundaryElementsCardArray)
# #         nExtBnd = 2*meshout.nElementBoundaries_global-sumCard
# #         nIntBnd = meshout.nElementBoundaries_global-nExtBnd
# #         meshout.nExteriorElementBoundaries_global=nExtBnd
# #         meshout.nInteriorElementBoundaries_global=nIntBnd

# #         #global edge numbers on exterior boundary
# #         meshout.exteriorElementBoundariesArray=numpy.zeros((nExtBnd,),'i')
# #         meshout.interiorElementBoundariesArray=numpy.zeros((nIntBnd,),'i')

# #         #enumerate interior and exterior
# #         interiorI = 0
# #         exteriorI = 0
# #         for ebN in range(meshout.nElementBoundaries_global):
# #             if elementBoundaryElementsCardArray[ebN] == 1:
# #                 meshout.exteriorElementBoundariesArray[exteriorI]=ebN
# #                 exteriorI+= 1
# #             else:
# #                 meshout.interiorElementBoundariesArray[interiorI]=ebN
# #                 interiorI+= 1
# #             #end if on card
# #         #end ebN

# #         #now track which nodes are on the boundary
# #         #maps edge --> its nodes
# #         #this is just the edgelist in triangle
# #         meshout.elementBoundaryNodesArray = numpy.zeros((nEdges_global,spaceDim),
# #                                                           'i')
# #         for edgeN in range(nEdges_global):
# #             meshout.elementBoundaryNodesArray[edgeN,:spaceDim]=edgeArray[edgeN,:spaceDim]-nbase
# #         #end edgeN
# #         #2d so edgeList is same as elementBoundaryNodesArray
# #         meshout.edgeNodesArray = numpy.zeros((nEdges_global,spaceDim),
# #                                                'i')
# #         for edgeN in range(nEdges_global):
# #             meshout.edgeNodesArray[edgeN,:spaceDim]=edgeArray[edgeN,:spaceDim]-nbase
# #         #end edgeN

# #         #now tell mesh that it doesn't have the list interface
# #         meshout.hasListInterface=False
# #         #compute diameters array manually
# #         meshout.elementDiametersArray=numpy.zeros((meshout.nElements_global,),
# #                                                     'd')
# #         import math
# #         for eN in range(meshout.nElements_global):
# #             elen = numpy.zeros((meshout.nElementBoundaries_element,),
# #                                  'd')
# #             for eloc in range(meshout.nElementBoundaries_element):
# #                 eg = meshout.elementBoundariesArray[eN,eloc] #glocal edge number
# #                 n0,n1 = meshout.elementBoundaryNodesArray[eg,0:2] #global node numbers number
# #                 de = meshout.nodeArray[n0,:]-meshout.nodeArray[n1,:]
# #                 elen[eloc] = math.sqrt(de[0]**2 + de[1]**2)
# #             #end eloc
# #             meshout.elementDiametersArray[eN]=max(elen)
# #         #end eN
# #         meshout.hasGeometricInfo=True
#         return meshout
    #end convertToProteusMesh
    def convertFromProteusMesh(self,meshin,verbose=0):
        """
        generate a Triangle mesh representation from an proteus mesh.
        This version will copy over the nodes and elements from the
        proteus mesh.

        Deletes the existing Triangle mesh and regenerates.
        """

        #first check that the input mesh has the correct dimensions and
        #minimal information necessary
        spaceDim = 2
        assert(meshin.nElementBoundaries_element == spaceDim+1)
        assert(not meshin.nodeArray == None)
        assert(not meshin.elementNodesArray == None)

        #get a clean slate
        tri0 = triangleWrappers.new()
        #don't set any markers for now
        #input array should be nNodes by spacedim
        nodesIn = meshin.nodeArray[:,:spaceDim]
        triangleWrappers.setPoints(tri0,nodesIn)
        #triangle array should be nElements x 3
        triAin  = meshin.elementNodesArray[:,:spaceDim+1]
        triangleWrappers.setTriangles(tri0,triAin)

        flagsAdd = "r" #use element and node connections
        flags = flagsAdd+self.baseFlags
        if flags.find('v') >= 0:
            self.makeVoronoi = True
            if verbose > 0:
                print 'readNodeAndEle makeVoronoi= ',self.makeVoronoi,' flags= ',flags
        #end if
        self.lastFlags = flags
        #reset data representations if necessary
        self.resetDefaultTrirep(index=0,verbose=verbose)
        if self.makeVoronoi:
            self.resetDefaultVorrep(index=0,verbose=verbose)
        #end vorrep initialized

        if self.makeVoronoi:
            triangleWrappers.applyTriangulate(flags,tri0,self.trirep[0],self.vorrep[0])
            #handles no longer contain trivial representations
            self.trirepDefaultInit[0] = False
            self.vorrepDefaultInit[0] = False

        else:
            triangleWrappers.applyTriangulateNoVoronoi(flags,tri0,self.trirep[0])
            #handle no longer contains trivial representations
            self.trirepDefaultInit[0] = False
        #end if

        #clean up explicitly
        del tri0
    #end convertFromProteusMesh
    ##################################################
    #read from file routines
    ##################################################
    def readFromNodeFile(self,filebase,flagsAdd="",verbose=0):
        """
        read triangle representation from filebase.node
        files. assumes the nbase data member is set appropriately
        """
        reader = TriangleUtils.TriangleInputFileReader()
        if verbose > 0:
            print 'readFromNodeAndEleFiles filebase= ',filebase
        #end if
        #could specify comment character too
        nodeDataInfo,nodeData = reader.readNodes(filebase,nbase=self.nbase)
        nodes,nodesA,nodesM = (nodeData['nodes'],
                               nodeData['nodeAttributes'],
                               nodeData['nodeMarkers'])

        if verbose > 3:
            print 'Nodes: nodeInfo= ',nodeDataInfo
            print """Nodes: nodes= \n%s\n nodeAttributes= \n%s\n nodeMarkers= \n%s
            """ % (nodes,nodesA,nodesM)
        #end if

        #now create an initial representation
        tri0 = triangleWrappers.new()

        if nodesM == None:
            triangleWrappers.setPoints(tri0,nodes)
        else:
            triangleWrappers.setPointsAndMarkers(tri0,nodes,nodesM)
        if not nodesA == None:
            triangleWrappers.setPointAttributes(tri0,nodesA)

        #run triangleWrappers on it using the base flags and whatever else was
        #added
        flags = flagsAdd+self.baseFlags
        if flags.find('v') >= 0:
            self.makeVoronoi = True
            if verbose > 0:
                print 'readNodeAndEle makeVoronoi= ',self.makeVoronoi,' flags= ',flags
        #end if
        self.lastFlags = flags

        #reset data representations if necessary
        self.resetDefaultTrirep(index=0,verbose=verbose)
        if self.makeVoronoi:
            self.resetDefaultVorrep(index=0,verbose=verbose)
        #end vorrep initialized


        if self.makeVoronoi:
            triangleWrappers.applyTriangulate(flags,tri0,self.trirep[0],self.vorrep[0])
            #handles no longer contain trivial representations
            self.trirepDefaultInit[0] = False
            self.vorrepDefaultInit[0] = False
        else:
            triangleWrappers.applyTriangulateNoVoronoi(flags,tri0,self.trirep[0])
            #handle no longer contains trivial representations
            self.trirepDefaultInit[0] = False
        #end if

        #do I need to clean up explicitly?
        del tri0
    def readFromNodeAndEleFiles(self,filebase,flagsAdd="",verbose=0):
        """
        read triangle representation from filebase.node and filebase.ele
        files. assumes the nbase data member is set appropriately
        """
        reader = TriangleUtils.TriangleInputFileReader()
        if verbose > 0:
            print 'readFromNodeAndEleFiles filebase= ',filebase
        #end if
        nodeDataInfo,nodeData = reader.readNodes(filebase,nbase=self.nbase)
        nodes,nodesA,nodesM = (nodeData['nodes'],
                               nodeData['nodeAttributes'],
                               nodeData['nodeMarkers'])

        if verbose > 3:
            print 'Nodes: nodeInfo= ',nodeDataInfo
            print """Nodes: nodes= \n%s\n nodeAttributes= \n%s\n nodeMarkers= \n%s
            """ % (nodes,nodesA,nodesM)
        #end if

        triDataInfo,triData = reader.readTriangles(filebase,nbase=self.nbase)
        triangles,trianglesA = triData['triangles'],triData['triangleAttributes']

        if verbose > 3:
            print 'Triangles: triInfo= ',triDataInfo
            print """Triangles: elems= \n%s\n triAttributes= \n%s
            """ % (triangles,trianglesA)
        #end if

        #now create an initial representation
        tri0 = triangleWrappers.new()

        if nodesM == None:
            triangleWrappers.setPoints(tri0,nodes)
        else:
            triangleWrappers.setPointsAndMarkers(tri0,nodes,nodesM)
        if not nodesA == None:
            triangleWrappers.setPointAttributes(tri0,nodesA)
        #end if

        triangleWrappers.setTriangles(tri0,triangles)
        if not trianglesA == None:
            triangleWrappers.setTriangleAttributes(tri0,trianglesA)

        #run triangulate on it using the base flags and whatever else was
        #added
        flags = flagsAdd+self.baseFlags
        if flags.find('v') >= 0:
            self.makeVoronoi = True
            if verbose > 0:
                print 'readNodeAndEle makeVoronoi= ',self.makeVoronoi,' flags= ',flags
        #end if
        self.lastFlags = flags

        #reset data representations if necessary
        self.resetDefaultTrirep(index=0,verbose=verbose)
        if self.makeVoronoi:
            self.resetDefaultVorrep(index=0,verbose=verbose)
        #end vorrep initialized

        if self.makeVoronoi:
            triangleWrappers.applyTriangulate(flags,tri0,self.trirep[0],self.vorrep[0])
            #handles no longer contain trivial representations
            self.trirepDefaultInit[0] = False
            self.vorrepDefaultInit[0] = False
        else:
            triangleWrappers.applyTriangulateNoVoronoi(flags,tri0,self.trirep[0])
            #handle no longer contains trivial representations
            self.trirepDefaultInit[0] = False
        #end if

        #do I need to clean up explicitly?
        del tri0
    def readFromPolyFile(self,filebase,flagsAdd="",verbose=0):
        """
        read triangle representation from filebase.poly file
        assumes the nbase data member is set appropriately
        """
        reader = TriangleUtils.TriangleInputFileReader()
        if verbose > 0:
            print 'readFromPolyFile filebase= ',filebase
        #end if

        polyDataInfo,polyData = reader.readPoly(filebase,nbase=self.nbase,
                                                verbose=verbose)

        if verbose > 3:
            for type in ['node','segment','hole','region']:
                print 'Poly: ',type,'Info= \n',polyDataInfo[type],'\n'
            #end for
        #end if

        nodes,nodesA,nodesM = (polyData['node']['nodes'],
                               polyData['node']['nodeAttributes'],
                               polyData['node']['nodeMarkers'])
        segments,segmentsM   = (polyData['segment']['segments'],
                                polyData['segment']['segmentMarkers'])
        holes                 = polyData['hole']['holes']
        regions               = polyData['region']['regions']

        if verbose > 3:
            print """Poly file read:
            nodes   = \n%s\n nodeAttributes= \n%s\n nodeMarkers= \n%s\n
            segments= \n%s\n segmentMarkers= \n%s\n
            holes   = \n%s\n
            regions = \n%s\n
            """ % (nodes,nodesA,nodesM,segments,segmentsM,holes,regions)
        #end if
        tri0 = triangleWrappers.new()

        if nodesM == None:
            triangleWrappers.setPoints(tri0,nodes)
        else:
            nodesMtmp = numpy.zeros(nodesM.shape,'i')
            nodesMtmp[:] = nodesM
            nodesM = nodesMtmp
            triangleWrappers.setPointsAndMarkers(tri0,nodes,nodesM)
        if not nodesA == None:
            triangleWrappers.setPointAttributes(tri0,nodesA)
        #end if
        if segmentsM == None:
            triangleWrappers.setSegments(tri0,segments)
        else:
            triangleWrappers.setSegmentsAndMarkers(tri0,segments,segmentsM)
        #end if
        if (not holes == None):
            triangleWrappers.setHoles(tri0,holes)
        #end if
        if (not regions == None):
            #print 'setting trin1 regions=\n',regions2
            triangleWrappers.setRegions(tri0,regions)
        #end if

        flags = flagsAdd+self.baseFlags
        if flags.find('p') == -1:
            print 'flags = ',flags,' must have p, appending'
            flags += "p"
        #end
        if flags.find('v') >= 0:
            self.makeVoronoi = True
            if verbose > 0:
                print 'readNodeAndEle makeVoronoi= ',self.makeVoronoi,' flags= ',flags
        #end if
        self.lastFlags = flags
        #reset data representations if necessary
        self.resetDefaultTrirep(index=0,verbose=verbose)
        if self.makeVoronoi:
            self.resetDefaultVorrep(index=0,verbose=verbose)
        #end vorrep initialized

        if self.makeVoronoi:
            triangleWrappers.applyTriangulate(flags,tri0,self.trirep[0],self.vorrep[0])
            #handles no longer contain trivial representations
            self.trirepDefaultInit[0] = False
            self.vorrepDefaultInit[0] = False
        else:
            triangleWrappers.applyTriangulateNoVoronoi(flags,tri0,self.trirep[0])
            #handle no longer contains trivial representations
            self.trirepDefaultInit[0] = False
        #end if

        #clean up?
        del tri0
    ##################################################
    #output routines
    ##################################################
    def writeToFile(self,filebase,verbose=0):
        """
        Just write out basic files for triangulateion
        Still need to write out Voronoi diagram
        """
        TriangleUtils.writeOutTriangulation(self.trirep[0],filebase,
                                            self.nbase,verbose)

    #end def

    def viewShowme(self):
        """
        just call showme for the mesh held in rep, uses meshshow-tmp file
        """
        filebase="meshshow-tmp"
        self.writeToFile(filebase)
        globshowme= TriangleUtils.showmeCmdBase
        showmecmd = """%s  %s """ % (globshowme,filebase)

        failure = 0
        failure = TriangleFileUtils.checkFileExists(globshowme)
        failure = os.system(showmecmd)

        return failure

    #end viewShowme
#end TriangleBaseMesh

########################################################################
#define some simple helper functions and examples on calling the code
########################################################################

def TriangleCall3(filebase="trimesh",baseFlags="zen",
                  flagsAdd="",verbose=0):
    import TriangleIface
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

    print 'viewing mesh generated'
    mesh.viewShowme()

#end TriangleCall3
def exProteusMesh0(filebase="trimesh",baseFlags="zen",
                     flagsAdd="",viewMesh=1,verbose=0):
    """
    create a Triangle mesh representation
    read it in from a file and initialize
    convert to an proteusMesh
    """
    import TriangleIface
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
    proteusMesh = mesh.convertToProteusMesh(verbose)
    proteusMesh.writeEdgesGnuplot2("gnuMesh") #uses array interface
    if viewMesh > 0:
        proteusMesh.viewMeshGnuplotPipe("gnuMesh")
    matmesh = "matlabMesh"
    proteusMesh.buildMatlabMeshDataStructures(matmesh)

#end exAdh
def exProteusMesh1(baseFlags="zen",viewMesh=1,verbose=0):
    """
    create an proteusMesh for a rectangle
    create a Triangle mesh representation
    convert the proteusMesh mesh to the Triangle mesh representation
    """
    import MeshTools
    import TriangleIface
    #simple domain for now
    Lx = 2.0
    Ly = 1.0
    nx = 11
    ny = 6
    proteusMesh = MeshTools.TriangularMesh()
    proteusMesh.constructTriangularMeshOnRectangle(Lx,Ly,nx,ny)
    proteusMesh.writeEdgesGnuplot2("gnuMesh") #uses array interface
    if viewMesh > 0:
        proteusMesh.viewMeshGnuplotPipe("gnuMesh")

    nbase = 0
    if baseFlags.find('z') == -1:
        nbase=1
    #end if
    trimesh = TriangleIface.TriangleBaseMesh(baseFlags=baseFlags,
                                             nbase=nbase,
                                             verbose=verbose)
    trimesh.convertFromProteusMesh(proteusMesh,verbose=verbose)

    if viewMesh > 0:
        print 'viewing mesh generated from proteus mesh'
        trimesh.viewShowme()
    #end if

#end exAdh

def exProteusLaplace1(filebase="trimesh",baseFlags="zen",
                        flagsAdd="",viewMesh=1,verbose=0):
    import QuadTools
    import FemTools
    import TriangleIface
    import PoissonTestProblems
    import ScalarTransport
    import TimeIntegrationTools
    import LinearAlgebraTools
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
    proteusMesh = mesh.convertToProteusMesh(verbose)
    proteusMesh.writeEdgesGnuplot2("gnuMesh") #uses array interface
    if viewMesh > 0:
        proteusMesh.viewMeshGnuplotPipe("gnuMesh")
    matmesh = "matlabMesh"
    proteusMesh.buildMatlabMeshDataStructures(matmesh)

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
    A  = numpy.zeros((2,2),'d')
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
            FemSpace = FemTools.C0_AffineLinearOnSimplexWithNodalBasis(proteusMesh,nd)
        else:
            FemSpace = FemTools.C0_AffineQuadraticOnSimplexWithNodalBasis(proteusMesh,nd)
    else:
        if order == 1:
            FemSpace = FemTools.DG_AffineLinearOnSimplexWithNodalBasis(proteusMesh,nd)
        else:
            FemSpace = FemTools.DG_AffineQuadraticOnSimplexWithNodalBasis(proteusMesh,nd)

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
    y  = numpy.zeros((system.dim,),'d')
    dy = numpy.zeros((system.dim,),'d')
    r  = numpy.zeros((system.dim,),'d')
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
    nlSolver.solve(y,r)


    if verbose > 0:
        print 'nonlinear solver done'
    #end if

    #print out solution
    FemSpace.writeFunctionMatlab(u,"ex1soln")
    #can view results in matlab with
    #matlabMesh; ex1soln; pdeplot(p,e,t,'xydata',u,'contour','on')

#end exAdh



if __name__ == '__main__':
    #make sure python has proteus in path
    import os,sys
    from optparse import OptionParser
    parser = OptionParser()

    #options controlling simulation behavior
    parser.add_option('-P','--proteusDir',
                      default='/Users/mfarthin/Public/code/proteus-trunk/src',
                      help="""where to find proteus library""")
    parser.add_option('--baseFlags',
                      default='zen',
                      help="""base flags for creating Triangle meshes""")
    parser.add_option('--filebase',
                      default='trimesh',
                      help="""base filename for reading a Triangle mesh""")
    parser.add_option('--tribase',
                      default='/Users/mfarthin/Public/code/proteus-packages/triangle/',
                      help="""base location for Triangle""")
    parser.add_option('--flagsAdd',
                      default='',
                      help="""flags to add when reading/creating Triangle meshes""")

    parser.add_option('-t','--testNum',
                      default=0,
                      help="""which test to use [0]
                      0 --- TriangleCall3 trimesh fro  .nodes and .ele file,
                            .poly file,  and .node file
                      1 --- TriangleCall3 with spiral
                      2 --- exProteusMesh0 with spiral
                      3 --- exProteusMesh0 with la
                      4 --- exProteusMesh0 with user specified file, and flags
                      5 --- exProteusMesh1 convert proteusMesh for rectangle
                      6 --- solve simple laplace equation with dirichlet bcs
                            on the la mesh
                      """)
    parser.add_option('-v','--verbose',
                      default=0,
                      help="""level of verbosity in simulation [0]""")

    #get options
    options, args = parser.parse_args(sys.argv[1:]) #get command line args

    #
    verbose  = int(options.verbose)
    testNum  = int(options.testNum)
    proteusDir = str(options.proteusDir)
    flagsAdd = str(options.flagsAdd)
    baseFlags= str(options.baseFlags)
    filebase = str(options.filebase)
    tribase  = str(options.tribase)

    sys.path.insert(0,proteusDir)

    if testNum == 0:
        #look at node and ele files
        filebase=os.path.join(tribase,"examples/trimesh")
        flagsAdd="r"
        TriangleCall3(filebase=filebase,flagsAdd=flagsAdd,
                      verbose=verbose)
        #generat from poly file
        filebase=os.path.join(tribase,"examples/trimesh")
        flagsAdd="p"
        TriangleCall3(filebase=filebase,flagsAdd=flagsAdd,
                      verbose=verbose)
        #generat from node file
        filebase=os.path.join(tribase,"examples/trimesh")
        flagsAdd=""
        TriangleCall3(filebase=filebase,flagsAdd=flagsAdd,
                      verbose=verbose)
    elif testNum == 1:
        filebase=os.path.join(tribase,"examples/spiral")
        baseFlags="en" #base 1?
        flagsAdd =""   #just node
        TriangleCall3(filebase=filebase,flagsAdd=flagsAdd,
                      baseFlags=baseFlags,verbose=verbose)
    elif testNum == 2:
        filebase=os.path.join(tribase,"examples/spiral")
        baseFlags="en" #base 1?
        flagsAdd =""   #just node
        exProteusMesh0(filebase=filebase,flagsAdd=flagsAdd,
                         baseFlags=baseFlags,verbose=verbose)
    elif testNum == 3:
        filebase=os.path.join(tribase,"examples/la")
        baseFlags="en" #base 1?
        flagsAdd ="pa"   #from poly file
        exProteusMesh0(filebase=filebase,flagsAdd=flagsAdd,
                         baseFlags=baseFlags,verbose=verbose)
    elif testNum == 4:
        exProteusMesh0(filebase=filebase,flagsAdd=flagsAdd,
                         baseFlags=baseFlags,verbose=verbose)
    elif testNum == 5:
        exProteusMesh1(baseFlags=baseFlags,verbose=verbose)
    elif testNum == 6:
        filebase=os.path.join(tribase,"examples/la")
        baseFlags="en" #base 1?
        flagsAdd ="pqa"   #from poly file
        exProteusLaplace1(filebase,baseFlags,
                            flagsAdd=flagsAdd,viewMesh=1,verbose=0)
        #can view results in matlab with
        #matlabMesh; ex1soln; pdeplot(p,e,t,'xydata',u,'contour','on')
    else:
        print 'PROBLEM testNum= ',testNum,' not recognized'
    #end if

def testGenerateTriangulationFromPointSet(points):
    """
    input: point-set

    generate triangle representation of the points

    output: an array contaning the points, and element to Node connectivity
    """
    #mwf debug
    #import pdb
    #pdb.set_trace()
    #default representation for holding input points
    tri0 = triangleWrappers.new()

    triangleWrappers.setPoints(tri0,points[:,:2])

    flags = "qz" #just use simple quality mesh generation, number from zero

    #default representation for holding output points
    tri1 = triangleWrappers.new()

    triangleWrappers.applyTriangulateNoVoronoi(flags,tri0,tri1)

    #if returning nodeArray and elementNodesArray outside of routine
    #use getPointsCopy, getTrianglesCopy
    #otherwise if doing all the necessary generation of quadrature points
    #and weights internally don't need deep copy
    nodeArray2d = triangleWrappers.getPointsCopy(tri1)
    elementNodesArray = triangleWrappers.getTrianglesCopy(tri1)



    #clean up
    del tri0
    del tri1

    return nodeArray2d,elementNodesArray

def testGenerateSSIPtriangulation(points):
    """
    input: vertices and SSIP points belonging to a single element

    generate triangle representation of the points

    output: an array contaning the points, and element quadrature weights for the points

    test with input
import numpy
from proteus import TriangleTools
points = numpy.array([[0.0,0.0,0.0],[0.5,0.4,0.],[1.0,0.0,0.0],[0.2,0.3,0.0],[0.0,1.0,0.0]])
dV,x = TriangleTools.testGenerateSSIPtriangulation(points)

    """
    #mwf debug
    #import pdb
    #pdb.set_trace()
    #default representation for holding input points
    tri0 = triangleWrappers.new()

    triangleWrappers.setPoints(tri0,points[:,:2])

    flags = "z"#"qz" #just use simple quality mesh generation, number from zero

    #default representation for holding output points
    tri1 = triangleWrappers.new()

    triangleWrappers.applyTriangulateNoVoronoi(flags,tri0,tri1)

    #doing all the necessary generation of quadrature points
    #and weights internally don't need deep copy
    nodeArray2d = triangleWrappers.getPoints(tri1)
    elementNodesArray = triangleWrappers.getTriangles(tri1)

    import Quadrature,cfemIntegrals
    #try Gauss quadrature different orders
    #subElementQuadrature = Quadrature.GaussTriangle()
    #subElementQuadrature.setOrder(2)  #order(1)
    #What about vertex quadrature
    subElementQuadrature = Quadrature.LobattoTriangle()
    subElementQuadrature.setOrder(3)

    nSubQuadraturePoints = len(subElementQuadrature.weights)
    nElements = elementNodesArray.shape[0]
    nQuadraturePoints = nElements*nSubQuadraturePoints
    nSpace = 2
    nDOF_element = 3
    weights_ref = numpy.array(subElementQuadrature.weights)
    points_ref  = numpy.array(subElementQuadrature.points)
    #loop through elements, compute area, map reference points to physical space,
    #and compute physical space weights

    #to be consistent need to have nodes be 3d
    nodeArray = numpy.zeros((nodeArray2d.shape[0],3),'d')
    nodeArray[:,0:2] = nodeArray2d
    #need to store psi and grad_psi to avoid recalculating across a lot of elements?
    #shape functions  are 1-x-y, x, y
    psi = numpy.zeros((nSubQuadraturePoints,nDOF_element),'d')
    psi[:,0] = 1.0-points_ref[:,0]-points_ref[:,1]
    psi[:,1] = points_ref[:,0]
    psi[:,2] = points_ref[:,1]
    #for k in range(nSubQuadraturePoints):
    #    psi[k,0] = 1.0-subElementQuadrature.points[k][0]-subElementQuadrature.points[k][1]
    #    psi[k,1] = subElementQuadrature.points[k][0]
    #    psi[k,2] = subElementQuadrature.points[k][1]
    grad_psi = numpy.zeros((nSubQuadraturePoints,nDOF_element,nSpace),'d')
    #shape functions  are 1-x-y, x, y
    grad_psi[:,0,0].fill(-1.); grad_psi[:,0,1].fill(-1.)
    grad_psi[:,1,0].fill(1.0);
    grad_psi[:,2,1].fill(1.0)
    #waste space to reuse code
    jacobianArray    = numpy.zeros((nElements,nSubQuadraturePoints,nSpace,nSpace),'d')
    jacobianDetArray = numpy.zeros((nElements,nSubQuadraturePoints),'d')
    jacobianInvArray = numpy.zeros((nElements,nSubQuadraturePoints,nSpace,nSpace),'d')

    cfemIntegrals.parametricMaps_getJacobianValues(grad_psi,
                                                   elementNodesArray,
                                                   nodeArray,
                                                   jacobianArray,
                                                   jacobianDetArray,
                                                   jacobianInvArray)


    jacobianDetArrayAbs = numpy.absolute(jacobianDetArray)
    #weights and points shaped like nSubElements x nQuadraturePointsSubElement
    q_sub_dV = numpy.zeros((nElements,nSubQuadraturePoints),'d')
    q_sub_x  = numpy.zeros((nElements,nSubQuadraturePoints,3),'d')

    cfemIntegrals.calculateIntegrationWeights(jacobianDetArrayAbs,
                                              weights_ref,
                                              q_sub_dV)

    cfemIntegrals.parametricMaps_getValues(psi,
                                           elementNodesArray,
                                           nodeArray,
                                           q_sub_x)


    #clean up
    del tri0
    del tri1

    #exit
    return q_sub_dV,q_sub_x
