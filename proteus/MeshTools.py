"""
Tools for creating and manipulating 1,2, and 3D meshes.

.. inheritance-diagram:: proteus.MeshTools
   :parts: 1
"""
from EGeometry import *
import numpy as np
import array
from Archiver import *
from LinearAlgebraTools import ParVec_petsc4py
from .Profiling import logEvent,memory

class Node:
    """A numbered point in 3D Euclidean space

    :ivar N: node number
    :ivar p: Euclidean coordinates

    Comparison operators and a hash value are defined using the 3-tuple of
    coordinates.  This allows using Node objects and tuples of node objects as
    dictionary keys, but in that use case one should be careful not to modify
    the node coordinates.

    >>> n0 = Node(nodeNumber=0,x=0.0,y=0.0,z=0.0)
    >>> n1 = Node(nodeNumber=1,x=1.0,y=1.0,z=1.0)
    >>> n1 >= n0
    True
    """
    xUnitVector = EVec(1.0,0.0,0.0)
    yUnitVector = EVec(0.0,1.0,0.0)
    zUnitVector = EVec(0.0,0.0,1.0)
    def __init__(self,nodeNumber=0,x=0.0,y=0.0,z=0.0):
        self.N=nodeNumber
        self.p=EVec(x,y,z)
        self.basis = [Node.xUnitVector,
                      Node.yUnitVector,
                      Node.zUnitVector]
        self.elementBoundaries=[]
        self.barycenter = self.p
        self.length = 1.0
        self.diameter=self.length
        self.innerDiameter=self.length
        self.hasGeometricInfo = True
        self.unitNormal = Node.xUnitVector
        self.nodes=(self,)
    def computeGeometricInfo(self):
        pass
    def __str__(self):
        return str(self.N)+":"+str(self.p)
    def __hash__(self):
        return hash((self.p[X],self.p[Y],self.p[Z]))
    def __lt__(self,other):
        return  (self.p[X],self.p[Y],self.p[Z]) < \
               (other.p[X],other.p[Y],other.p[Z])
    def __le__(self,other):
        return  (self.p[X],self.p[Y],self.p[Z]) <= \
               (other.p[X],other.p[Y],other.p[Z])
    def __eq__(self,other):
        return  (self.p[X],self.p[Y],self.p[Z]) == \
               (other.p[X],other.p[Y],other.p[Z])
    def __ne__(self,other):
        return  (self.p[X],self.p[Y],self.p[Z]) != \
               (other.p[X],other.p[Y],other.p[Z])
    def __gt__(self,other):
        return  (self.p[X],self.p[Y],self.p[Z]) > \
               (other.p[X],other.p[Y],other.p[Z])
    def __ge__(self,other):
        return  (self.p[X],self.p[Y],self.p[Z]) >= \
               (other.p[X],other.p[Y],other.p[Z])

class Element:
    """An numbered polytope in R^n

    :ivar N: element number
    :ivar nodes: sorted tuple of nodes defining the polytope

    The nodes data member can be used as a dictionary key for the polytope as
    long as the nodes aren't later modified.
    """
    def __init__(self,elementNumber=0,nodes=[]):
        self.N = elementNumber
        nodeList = nodes
        nodeList.sort()
        self.nodes = tuple(nodeList)
        self.elementBoundaries=[]

class Edge(Element):
    xUnitVector = EVec(1.0,1.0,0.0)
    yUnitVector = EVec(0.0,1.0,0.0)
    zUnitVector = EVec(0.0,0.0,1.0)
    """
    1D Element--a line connecting two Nodes

    The nodes are stored as a lexicographically sorted node list.
    """
    def __init__(self,edgeNumber=0,nodes=[]):
        #Element.__init__(self,edgeNumber,nodes)
        #inline Element.__init__
        self.N = edgeNumber
        nodeList = nodes
        nodeList.sort()
        self.nodes = tuple(nodeList)
        #self.nodes=nodes
        #self.nodes=nodes[:]
        #self.nodes.sort()
        self.elementBoundaries = [self.nodes[1],self.nodes[0]]
        self.hasGeometricInfo = False

    def computeGeometricInfo(self):
        if not self.hasGeometricInfo:
            self.basis = [self.nodes[1].p - self.nodes[0].p,
                          Edge.yUnitVector,
                          Edge.zUnitVector]
            self.barycenter = (self.nodes[0].p + self.nodes[1].p)/2.0
            self.length = enorm(self.basis[0])
            self.normal = EVec(-self.basis[0][Y], self.basis[0][X],0.0)
            norm = enorm(self.normal)
            if  norm:
                self.unitNormal = self.normal/norm
            else:
                #in 3D edge normals don't make sense in general so above
                #may divide by zero if edge has zero projection onto x-y plane
                self.normal = EVec(0.0, -self.basis[0][Z], self.basis[0][Y])
                self.unitNormal = self.normal/enorm(self.normal)
            self.diameter=self.length
            self.innerDiameter = self.length
            self.hasGeometricInfo = True
            self.nodeUnitNormalList=[]
            self.nodeUnitNormalList.append(-self.basis[0]/self.length)
            self.nodeUnitNormalList.append(self.basis[0]/self.length)
            self.elementBoundaryUnitNormalList=self.nodeUnitNormalList
            self.elementBoundaryJacobianList=[Edge.xUnitVector,Edge.xUnitVector]
def getNodesFromEdges(edges):
    """Extract the subset of nodes from a list of edges."""
    nodes={}
    for e in edges:
        for n in e.nodes:
            nodes[n]=n
    return nodes.values()

class Polygon(Element):
    """An abstract 2D element--a closed set of Edges connecting a set of Nodes.

    The nodes and edges are stored as lexicographically sorted lists."""
    def __init__(self,polygonNumber=0,nodes=[]):
        Element.__init__(self,polygonNumber,nodes)
        #the edges have to be set up by the specific polygon
        self.edges=[]

def getEdgesFromPolygons(polygons):
    """Extract the subset of edges from a list of polygons"""
    edges={}
    for p in polygons:
        for e in p.edges:
            edges[e.nodes] = e
    return edges.values()

class Triangle(Polygon):
    """A 2D triangular element"""
    edgeMap = {(1,2):0,(0,2):1,(0,1):2}
    zUnitVector = EVec(0.0,0.0,1.0)
    def __init__(self,triangleNumber=0,nodes=[],edgeDict=None):
        #Polygon.__init__(self,triangleNumber,nodes)
        #inline
        self.edges=[]
        #Element.__init__
        #inline
        self.N = triangleNumber
        nodeList = nodes
        nodeList.sort()
        self.nodes = tuple(nodeList)
        #self.nodes=nodes[:]
        #self.nodes.sort()
        #
        edgeNodeList = [(self.nodes[1],self.nodes[2]),
                        (self.nodes[0],self.nodes[2]),
                        (self.nodes[0],self.nodes[1])]
        if edgeDict is None:
            self.edges = [Edge(eN,list(edgeNodes)) for \
                          eN,edgeNodes in enumerate(edgeNodeList)]
        else:
            self.edges = [edgeDict[edgeNodes] for edgeNodes in edgeNodeList]
        self.hasGeometricInfo=False
        self.elementBoundaries=self.edges

    def computeGeometricInfo(self):
        if not self.hasGeometricInfo:
            self.barycenter = (self.nodes[0].p +
                               self.nodes[1].p +
                               self.nodes[2].p)/3.0
            self.basis = [ n.p - self.nodes[0].p for n in self.nodes[1:]]
            self.basis.append(Triangle.zUnitVector)
            self.linearMap = ETen(self.basis[0],self.basis[1],self.basis[2])
            self.normal = ecross(self.basis[0],self.basis[1])
            normNormal = enorm(self.normal)
            self.unitNormal = self.normal/normNormal
            self.area = 0.5*normNormal
            for e in self.edges: e.computeGeometricInfo()
            self.diameter = max([e.length for e in self.edges])
            self.innerDiameter = 4.0*self.area/sum(
                [e.length for e in self.edges])
            self.edgeUnitNormalList=[]
            for nNt,eN in Triangle.edgeMap.iteritems():
                unitNormal = self.edges[eN].unitNormal
                if edot(unitNormal,self.nodes[nNt[0]].p - self.nodes[eN].p) < 0:
                    unitNormal*=-1.0
                self.edgeUnitNormalList.append(unitNormal)
            self.elementBoundaryUnitNormalList = self.edgeUnitNormalList
            self.hasGeometricInfo=True

class Quadrilateral(Polygon):
    """A 2D quadrilateral element"""
    def __init__(self,quadrilateralNumber=0,edges=[],simple=True):
        Polygon.__init__(self,quadrilateralNumber)
        self.edges = edges
        nodeList = getNodesFromEdges(self.edges)
        nodeList = self.sortNodes(nodeList)
        self.nodes = tuple(nodeList)
        self.hasGeometricInfo = False
        self.elementBoundaries = self.edges
        # This boolean flags whether the quadrilateral is simple
        # (eg. a rectangle).  Certain features are more difficult
        # to implement if this is not the case.
        self.simple = True

    def sortNodes(self,nodeList):
        newList = [None] * 4
        coordinate_list = [1,1,1]

        # initialize coordinate mins and maxs
        xMin = nodeList[0].p[X]
        xMax = nodeList[0].p[X]
        yMin = nodeList[0].p[Y]
        yMax = nodeList[0].p[Y]
        zMin = nodeList[0].p[Z]
        zMax = nodeList[0].p[Z]
        for node in nodeList:
            if xMin > node.p[X]:
                xMin = node.p[X]
            if xMax < node.p[X]:
                xMax = node.p[X]
            if yMin > node.p[Y]:
                yMin = node.p[Y]
            if yMax < node.p[Y]:
                yMax = node.p[Y]
            if zMin > node.p[Z]:
                zMin = node.p[Z]
            if zMax < node.p[Z]:
                zMax = node.p[Z]

        # indentify degenerate coordinate space.
        # NOTE - this is not entirely accurate, but assumes
        # 2D quadrilateral objects are orthogonal to one of
        # the cononical coordinate axes

        if xMin==xMax:
            coordinate_list[0] = 0
        if yMin==yMax:
            coordinate_list[1] = 0
        if zMin==zMax:
            coordinate_list[2] = 0
        if sum(coordinate_list) !=2:
            assert 0, 'Invalid 2D quadrilateral object'

        for i, t in enumerate(coordinate_list):
            if t == 0:
                case = i

        # x is degenerate variable
        if case == 0:
            var1 = 1        # y marked as first node
            var2 = 2        # z marked as second
            var1_min = yMin
            var1_max = yMax
            var2_min = zMin
            var2_max = zMax
        # y is degenerate variable
        elif case == 1:
            var1 = 0        # x marked as first node
            var2 = 2        # z marked as second
            var1_min = xMin
            var1_max = xMax
            var2_min = zMin
            var2_max = zMax
        # z is degenerate variable
        elif case == 2:
            var1 = 0        # x marked as first node
            var2 = 1        # y marked as second
            var1_min = xMin
            var1_max = xMax
            var2_min = yMin
            var2_max = yMax
        else:
            assert 0, 'Invalide Quadrilateral Mesh Case'

        for node in nodeList:
            if node.p[var1]==var1_min and node.p[var2]==var2_min:
                newList[0] = node
            elif node.p[var1]==var1_min and node.p[var2]==var2_max:
                newList[1] = node
            elif node.p[var1]==var1_max and node.p[var2]==var2_max:
                newList[2] = node
            elif node.p[var1]==var1_max and node.p[var2]==var2_min:
                newList[3] = node

#        import pdb
#        pdb.set_trace()

        for i,item in enumerate(newList):
            if not newList[i]:
                assert 0,'Quadrialteral Mesh Generation Error '+`newList`+" i = "+`i`
        return newList

    def computeGeometricInfo(self):
        if not self.hasGeometricInfo:
            for e in self.edges: e.computeGeometricInfo()
            #the nodes must lie in a plane
            #use triangles to compute area
            #grab one triangle
            t0 = Triangle(0,list(self.nodes[0:3]))
            t0.computeGeometricInfo()
            #find the nodes  that lie on the new edge,diagonal0
            for et in t0.edges:
                edgeIsNew=True
                for e in self.edges:
                    if e.nodes == et.nodes:
                        edgeIsNew=False
                if edgeIsNew:
                    break
            diagonal0=et
            t1 = Triangle(0,[self.nodes[3],
                             diagonal0.nodes[0],
                             diagonal0.nodes[1]])
            t1.computeGeometricInfo()
            #get normal from one of the triangles
            self.unitNormal = t0.unitNormal
            self.area = t0.area + t1.area
            #find the long diagonal
            diagonalNode=0
            for n in self.nodes[0:3]:
                if n != diagonal0.nodes[0] and n != diagonal0.nodes[1]:
                    diagonalNode=n
                    break;
            diagonal1 = Edge(0,[n,self.nodes[3]])
            diagonal1.computeGeometricInfo()
            self.diameter = max(diagonal1.length,diagonal0.length)
            self.innerDiameter = 4.0*self.area/sum(
                [e.length for e in self.edges])
            # Calculate the coordinate of a simple quad
            if self.simple==True:
                self.xmin = self.nodes[0].p[X]
                self.ymin = self.nodes[0].p[Y]
                self.xmax = self.nodes[0].p[X]
                self.ymax = self.nodes[0].p[Y]
                for node in self.nodes:
                    if node.p[X] < self.xmin:
                        self.xmin = node.p[X]
                    elif node.p[X] > self.xmax:
                        self.xmax = node.p[X]
                    else:
                        pass
                    if node.p[Y] < self.ymin:
                        self.ymin = node.p[Y]
                    elif node.p[Y] > self.ymax:
                        self.ymax = node.p[Y]
                    else:
                        pass
                self.xmid = (self.xmin+self.xmax)/2.
                self.ymid = (self.ymin+self.ymax)/2.
                self.zmid = 0.

class Polyhedron(Element):
    """
    An abstract 3D Element--a closed set of Polygons connecting a set
    of Edges.

    The nodes and edges are stored as lexicographically sorted lists.
    """
    def __init__(self,polyhedronNumber=0,nodes=[]):
        Element.__init__(self,polyhedronNumber,nodes)
        self.edges=[]
        self.polygons=[]
    def __cmp__(self,other):
        return compareNodes(self.nodes,other.nodes)

class Tetrahedron(Polyhedron):
    """A 3D tetrahedral element"""
    triangleMap = {(1,2,3):0,(0,2,3):1,(0,1,3):2,(0,1,2):3}
    edgeMap = {(0,1): 0,
               (0,2): 1,
               (0,3): 2,
               (1,2): 3,
               (1,3): 4,
               (2,3): 5}

    def __init__(self,tetrahedronNumber,nodes,edgeDict=None,triangleDict=None):
        #Polyhedron.__init__(self,tetrahedronNumber,nodes)
        #inline
        #Element.__init__
        #inline
        self.N = tetrahedronNumber
        nodeList = nodes
        nodeList.sort()
        self.nodes = tuple(nodeList)
        #self.nodes=nodes[:]
        #self.nodes.sort()
        #
        triangleNodeList = [(self.nodes[1],
                             self.nodes[2],
                             self.nodes[3]),
                            (self.nodes[0],
                             self.nodes[2],
                             self.nodes[3]),
                            (self.nodes[0],
                             self.nodes[1],
                             self.nodes[3]),
                            (self.nodes[0],
                             self.nodes[1],
                             self.nodes[2])]
        if triangleDict is None:
            self.triangles = [Triangle(triangleNumber=tN,
                                       nodes=list(triangleNodes))
                              for tN,triangleNodes in
                              enumerate(triangleNodeList)]
        else:
            self.triangles = [triangleDict[triangleNodes] for triangleNodes in
                              triangleNodeList]
        self.polygons=self.triangles
        edgeNodeList = [(self.nodes[0],self.nodes[1]),
                        (self.nodes[0],self.nodes[2]),
                        (self.nodes[0],self.nodes[3]),
                        (self.nodes[1],self.nodes[2]),
                        (self.nodes[1],self.nodes[3]),
                        (self.nodes[2],self.nodes[3])]
        if edgeDict is None:
            self.edges = [Edge(edgeNumber=eN,nodes=list(edgeNodes)) for
                          eN,edgeNodes in enumerate(edgeNodeList)]
        else:
            self.edges = [edgeDict[edgeNodes] for edgeNodes in edgeNodeList]
        self.hasGeometricInfo=False
        self.elementBoundaries = self.triangles

    def computeGeometricInfo(self):
        if not self.hasGeometricInfo:
            for t in self.triangles: t.computeGeometricInfo()
            self.barycenter =(self.nodes[0].p +
                              self.nodes[1].p +
                              self.nodes[2].p +
                              self.nodes[3].p)/4.0
            self.basis = [n.p - self.nodes[0].p for n in self.nodes[1:]]
            self.linearMap = ETen(self.basis[0],self.basis[1],self.basis[2])
            self.volume = abs(edet(self.linearMap))/6.0
            self.diameter = max([t.diameter for t in self.triangles])
             #Zhang's formula for rho=innerDiameter of a simplex
            self.innerDiameter = 6.0*self.volume/sum([t.area for t in
                                                      self.triangles])
            self.triangleUnitNormalList=[]
            for nNt,tN in Tetrahedron.triangleMap.iteritems():
                unitNormal = self.triangles[tN].unitNormal
                if edot(unitNormal,self.nodes[nNt[0]].p - self.nodes[tN].p) < 0:
                    unitNormal *= -1.0
                self.triangleUnitNormalList.append(unitNormal)
            self.elementBoundaryUnitNormalList = self.triangleUnitNormalList
            self.hasGeometricInfo=True

class Hexahedron(Polyhedron):
    """A 3D hexahedral element"""
    def __init__(self,HN,quadrilaterals):
        Polyhedron.__init__(self,HN)
        self.N = HN
        self.quadrilaterals = quadrilaterals
        self.polygons = self.quadrilaterals
        self.edges = getEdgesFromPolygons(quadrilaterals)
        #self.nodes = getNodesFromEdges(self.edges)
        #self.nodes.sort()
        nodeList = getNodesFromEdges(self.edges)
        nodeList.sort()
        self.nodes = tuple(nodeList)
        self.hasGeometricInfo=False
        self.elementBoundaries = self.quadrilaterals

#todo add enum34 and replace with real Python enum
class MeshParallelPartitioningTypes:
    """
    fake an enum for parallel partitioning options
    """
    element = 0 ;  node     = 1

class Mesh:
    """A partition of a domain in R^n into elements.

    This is the base class for meshes. Contains routines for
    plotting the edges of the mesh in Matlab

    Attributes
    ----------
    elementBoundariesArray : array type
        This array lists the global edge number associated with every
        edge or face of an element.
    """
    #cek adding parallel support
    def __init__(self):
        #array interface
        self.nSubdomains_global=1
        self.sN = 0
        #node coordinates indexed by node number
        self.nNodes_global=0
        self.nNodes_subdomain=0
        self.nodeArray=None
        self.nodeVelocityArray=None
        self.nNodes_element=0
        #element node numbers, indexed by element number
        self.nElements_global=0
        self.nElements_proc=0
        self.elementNodesArray=None
        self.max_nElements_node=0
        self.nElements_node=None #mwf warning not calculated in buildPythonFromC
        self.nodeElementsArray=None
        self.nodeElementOffsets=None
        #element boundary numbers, indexed by element number
        self.nElementBoundaries_element=0
        self.elementBoundariesArray=None
        #element numbers, indexed by element boundary number and left(0) and right(1) element
        self.nElementBoundaries_global=0
        self.elementBoundaryElementsArray=None
        #local element boundary numbers, indexed by element boundary number and left(0) and right(1) element
        self.elementBoundaryLocalElementBoundariesArray=None
        #neighboring element  numbers, indexed by local element boundary number
        self.elementNeighborsArray=None
        #node numbers, indexed by element boundary number
        self.elementBoundaryNodesArray=None
        #element boundary numbers, indexed by interior/exterior
        #element boundary number
        self.interiorElementBoundariesArray=None
        self.nInteriorElementBoundaries_global=0
        self.exteriorElementBoundariesArray=None
        self.nExteriorElementBoundaries_global=0
        #edge node numbers, indexed by edge number
        self.nEdges_global=0
        self.edgeNodesArray=None
        self.nodeStarArray=None
        self.nodeStarOffsets=None
        self.h=0.0
        self.hMin=0.0
        self.hasGeometricInfo=False
        self.boundaryMesh=None
        #physical coordinates of element barycenters and elementBoundary barycenters
        self.elementBarycentersArray=None
        self.elementBoundaryBarycentersArray=None
        self.nodeDiametersArray=None
        self.nodeSupportArray=None
        #unique labels for classes of elements, elementBoundaries, nodes,
        self.elementMaterialTypes=None
        self.elementBoundaryMaterialTypes=None
        self.nodeMaterialTypes=None
        #parallel stuff
        self.nElements_owned=self.nElements_global
        self.nNodes_owned=self.nNodes_global
        self.nElementBoundaries_owned=self.nElementBoundaries_global
        self.nEdges_owned=self.nEdges_global
        self.elementOffsets_subdomain_owned=[0,self.nElements_global]
        self.elementNumbering_subdomain2global=np.arange(self.nElements_global,dtype='i')
        self.nodeOffsets_subdomain_owned=[0,self.nNodes_global]
        self.nodeNumbering_subdomain2global=np.arange(self.nNodes_global,dtype='i')
        self.elementBoundaryOffsets_subdomain_owned=[0,self.nElementBoundaries_global]
        self.elementBoundaryNumbering_subdomain2global=np.arange(self.nElementBoundaries_global,dtype='i')
        self.edgeOffsets_subdomain_owned=[0,self.nEdges_global]
        self.edgeNumbering_subdomain2global=np.arange(self.nEdges_global,dtype='i')
        self.subdomainMesh=self
        self.globalMesh = None
        self.arGridCollection=None
        self.arGrid=None
        self.nLayersOfOverlap = None
        self.parallelPartitioningType = MeshParallelPartitioningTypes.element
    def partitionMesh(self,nLayersOfOverlap=1,parallelPartitioningType=MeshParallelPartitioningTypes.node):
        import cmeshTools
        import Comm
        import flcbdfWrappers
        comm = Comm.get()
        self.comm=comm
        logEvent(memory("partitionMesh 1","MeshTools"),level=4)
        logEvent("Partitioning mesh among %d processors using partitioningType = %d" % (comm.size(),parallelPartitioningType))
        self.subdomainMesh=self.__class__()
        self.subdomainMesh.globalMesh = self
        self.subdomainMesh.cmesh=cmeshTools.CMesh()
        self.nLayersOfOverlap = nLayersOfOverlap; self.parallelPartitioningType = parallelPartitioningType
        logEvent(memory("partitionMesh 2","MeshTools"),level=4)
        if parallelPartitioningType == MeshParallelPartitioningTypes.node:
            #mwf for now always gives 1 layer of overlap
            (self.elementOffsets_subdomain_owned,
             self.elementNumbering_subdomain2global,
             self.nodeOffsets_subdomain_owned,
             self.nodeNumbering_subdomain2global,
             self.elementBoundaryOffsets_subdomain_owned,
             self.elementBoundaryNumbering_subdomain2global,
             self.edgeOffsets_subdomain_owned,
             self.edgeNumbering_subdomain2global) = flcbdfWrappers.partitionNodes(nLayersOfOverlap,self.cmesh,self.subdomainMesh.cmesh)
        else:
            (self.elementOffsets_subdomain_owned,
             self.elementNumbering_subdomain2global,
             self.nodeOffsets_subdomain_owned,
             self.nodeNumbering_subdomain2global,
             self.elementBoundaryOffsets_subdomain_owned,
             self.elementBoundaryNumbering_subdomain2global,
             self.edgeOffsets_subdomain_owned,
             self.edgeNumbering_subdomain2global) = flcbdfWrappers.partitionElements(nLayersOfOverlap,self.cmesh,self.subdomainMesh.cmesh)
        #
        logEvent(memory("partitionMesh 3","MeshTools"),level=4)
        self.subdomainMesh.buildFromC(self.subdomainMesh.cmesh)
        self.subdomainMesh.nElements_owned = self.elementOffsets_subdomain_owned[comm.rank()+1] - self.elementOffsets_subdomain_owned[comm.rank()]
        self.subdomainMesh.nNodes_owned = self.nodeOffsets_subdomain_owned[comm.rank()+1] - self.nodeOffsets_subdomain_owned[comm.rank()]
        self.subdomainMesh.nElementBoundaries_owned = self.elementBoundaryOffsets_subdomain_owned[comm.rank()+1] - self.elementBoundaryOffsets_subdomain_owned[comm.rank()]
        self.subdomainMesh.nEdges_owned = self.edgeOffsets_subdomain_owned[comm.rank()+1] - self.edgeOffsets_subdomain_owned[comm.rank()]

        comm.barrier()
        logEvent(memory("partitionMesh 4","MeshTools"),level=4)
        logEvent("Number of Subdomain Elements Owned= "+str(self.subdomainMesh.nElements_owned))
        logEvent("Number of Subdomain Elements = "+str(self.subdomainMesh.nElements_global))
        logEvent("Number of Subdomain Nodes Owned= "+str(self.subdomainMesh.nNodes_owned))
        logEvent("Number of Subdomain Nodes = "+str(self.subdomainMesh.nNodes_global))
        logEvent("Number of Subdomain elementBoundaries Owned= "+str(self.subdomainMesh.nElementBoundaries_owned))
        logEvent("Number of Subdomain elementBoundaries = "+str(self.subdomainMesh.nElementBoundaries_global))
        logEvent("Number of Subdomain Edges Owned= "+str(self.subdomainMesh.nEdges_owned))
        logEvent("Number of Subdomain Edges = "+str(self.subdomainMesh.nEdges_global))
        comm.barrier()
        logEvent("Finished partitioning")
        par_nodeDiametersArray = ParVec_petsc4py(self.subdomainMesh.nodeDiametersArray,
                                                 bs=1,
                                                 n=self.subdomainMesh.nNodes_owned,
                                                 N=self.nNodes_global,
                                                 nghosts=self.subdomainMesh.nNodes_global - self.subdomainMesh.nNodes_owned,
                                                 subdomain2global=self.nodeNumbering_subdomain2global)
        par_nodeDiametersArray.scatter_forward_insert()
        # comm.beginSequential()
        # from Profiling import memory
        # memory()
        # logEvent(memory("Partitioning Mesh","Mesh"),level=1)
        # del self.cmesh
        # #cmeshTools.deleteMeshDataStructures(self.cmesh)
        # logEvent(memory("Without global mesh","Mesh"),level=1)
        # comm.endSequential()
    def partitionMeshFromFiles(self,filebase,base,nLayersOfOverlap=1,parallelPartitioningType=MeshParallelPartitioningTypes.node):
        import cmeshTools
        import Comm
        import flcbdfWrappers
        comm = Comm.get()
        self.comm=comm
        logEvent(memory("partitionMesh 1","MeshTools"),level=4)
        logEvent("Partitioning mesh among %d processors using partitioningType = %d" % (comm.size(),parallelPartitioningType))
        self.subdomainMesh=self.__class__()
        self.subdomainMesh.globalMesh = self
        self.subdomainMesh.cmesh=cmeshTools.CMesh()
        self.nLayersOfOverlap = nLayersOfOverlap; self.parallelPartitioningType = parallelPartitioningType
        logEvent(memory("partitionMesh 2","MeshTools"),level=4)
        if parallelPartitioningType == MeshParallelPartitioningTypes.node:
            #mwf for now always gives 1 layer of overlap
            (self.elementOffsets_subdomain_owned,
             self.elementNumbering_subdomain2global,
             self.nodeOffsets_subdomain_owned,
             self.nodeNumbering_subdomain2global,
             self.elementBoundaryOffsets_subdomain_owned,
             self.elementBoundaryNumbering_subdomain2global,
             self.edgeOffsets_subdomain_owned,
             self.edgeNumbering_subdomain2global) = flcbdfWrappers.partitionNodesFromTetgenFiles(filebase,base,nLayersOfOverlap,self.cmesh,self.subdomainMesh.cmesh)
        else:
            (self.elementOffsets_subdomain_owned,
             self.elementNumbering_subdomain2global,
             self.nodeOffsets_subdomain_owned,
             self.nodeNumbering_subdomain2global,
             self.elementBoundaryOffsets_subdomain_owned,
             self.elementBoundaryNumbering_subdomain2global,
             self.edgeOffsets_subdomain_owned,
             self.edgeNumbering_subdomain2global) = flcbdfWrappers.partitionElementsFromTetgenFiles(filebase,base,nLayersOfOverlap,self.cmesh,self.subdomainMesh.cmesh)
        #
        logEvent(memory("partitionMesh 3","MeshTools"),level=4)
        self.buildFromCNoArrays(self.cmesh)
        self.subdomainMesh.buildFromC(self.subdomainMesh.cmesh)
        self.subdomainMesh.nElements_owned = self.elementOffsets_subdomain_owned[comm.rank()+1] - self.elementOffsets_subdomain_owned[comm.rank()]
        self.subdomainMesh.nNodes_owned = self.nodeOffsets_subdomain_owned[comm.rank()+1] - self.nodeOffsets_subdomain_owned[comm.rank()]
        self.subdomainMesh.nElementBoundaries_owned = self.elementBoundaryOffsets_subdomain_owned[comm.rank()+1] - self.elementBoundaryOffsets_subdomain_owned[comm.rank()]
        self.subdomainMesh.nEdges_owned = self.edgeOffsets_subdomain_owned[comm.rank()+1] - self.edgeOffsets_subdomain_owned[comm.rank()]

        comm.barrier()
        logEvent(memory("partitionMesh 4","MeshTools"),level=4)
        logEvent("Number of Subdomain Elements Owned= "+str(self.subdomainMesh.nElements_owned))
        logEvent("Number of Subdomain Elements = "+str(self.subdomainMesh.nElements_global))
        logEvent("Number of Subdomain Nodes Owned= "+str(self.subdomainMesh.nNodes_owned))
        logEvent("Number of Subdomain Nodes = "+str(self.subdomainMesh.nNodes_global))
        logEvent("Number of Subdomain elementBoundaries Owned= "+str(self.subdomainMesh.nElementBoundaries_owned))
        logEvent("Number of Subdomain elementBoundaries = "+str(self.subdomainMesh.nElementBoundaries_global))
        logEvent("Number of Subdomain Edges Owned= "+str(self.subdomainMesh.nEdges_owned))
        logEvent("Number of Subdomain Edges = "+str(self.subdomainMesh.nEdges_global))
        comm.barrier()
        logEvent("Finished partitioning")
        par_nodeDiametersArray = ParVec_petsc4py(self.subdomainMesh.nodeDiametersArray,
                                                 bs=1,
                                                 n=self.subdomainMesh.nNodes_owned,
                                                 N=self.nNodes_global,
                                                 nghosts=self.subdomainMesh.nNodes_global - self.subdomainMesh.nNodes_owned,
                                                 subdomain2global=self.nodeNumbering_subdomain2global)
        par_nodeDiametersArray.scatter_forward_insert()
        # comm.beginSequential()
        # from Profiling import memory
        # memory()
        # logEvent(memory("Partitioning Mesh","Mesh"),level=1)
        # del self.cmesh
        # #cmeshTools.deleteMeshDataStructures(self.cmesh)
        # logEvent(memory("Without global mesh","Mesh"),level=1)
        # comm.endSequential()
    def writeMeshXdmf(self,ar,name='',t=0.0,init=False,meshChanged=False,Xdmf_ElementTopology="Triangle",tCount=0, EB=False):
        if self.arGridCollection is not None:
            init = False
        elif not init:
            grids = ar.domain.findall("Grid")
            self.arGridCollection = grids[0]
            if EB:
                assert(len(grids) > 1)
                self.arEBGridCollection = grids[1]
        if init:
            self.arGridCollection = SubElement(ar.domain,"Grid",{"Name":"Mesh "+name,
                                                               "GridType":"Collection",
                                                               "CollectionType":"Temporal"})
            if EB:
                self.arEBGridCollection = SubElement(ar.domain,"Grid",{"Name":"EBMesh "+name,
                                                                       "GridType":"Collection",
                                                                       "CollectionType":"Temporal"})
        if self.arGrid is None or self.arTime.get('Value') != str(t):
            #
            #topology and geometry
            #
            if ar.global_sync:
                self.arGrid = SubElement(self.arGridCollection,"Grid",{"GridType":"Uniform"})
                self.arTime = SubElement(self.arGrid,"Time",{"Value":str(t),"Name":str(tCount)})
                topology = SubElement(self.arGrid,"Topology",
                                      {"Type":Xdmf_ElementTopology,
                                       "NumberOfElements":str(self.globalMesh.nElements_global)})
                elements = SubElement(topology,"DataItem",
                                      {"Format":ar.dataItemFormat,
                                       "DataType":"Int",
                                       "Dimensions":"%i %i" % (self.globalMesh.nElements_global,self.nNodes_element)})
                geometry = SubElement(self.arGrid,"Geometry",{"Type":"XYZ"})
                nodes    = SubElement(geometry,"DataItem",
                                      {"Format":ar.dataItemFormat,
                                       "DataType":"Float",
                                       "Precision":"8",
                                       "Dimensions":"%i %i" % (self.globalMesh.nNodes_global,3)})
                if ar.hdfFile is not None:
                    if ar.has_h5py:
                        elements.text = ar.hdfFilename+":/elements"+name+`tCount`
                        nodes.text = ar.hdfFilename+":/nodes"+name+`tCount`
                        if init or meshChanged:
                            ar.create_dataset_sync('elements'+name+`tCount`,
                                                    offsets=self.globalMesh.elementOffsets_subdomain_owned,
                                                    data=self.globalMesh.nodeNumbering_subdomain2global[self.elementNodesArray[:self.nElements_owned]])
                            ar.create_dataset_sync('nodes'+name+`tCount`,
                                                   offsets=self.globalMesh.nodeOffsets_subdomain_owned,
                                                   data=self.nodeArray[:self.nNodes_owned])
                    else:
                        assert False, "global_sync not supported  with pytables"
                else:
                    assert False, "global_sync not  supported with text heavy data"
            else:
                self.arGrid = SubElement(self.arGridCollection,"Grid",{"GridType":"Uniform"})
                self.arTime = SubElement(self.arGrid,"Time",{"Value":str(t),"Name":str(tCount)})
                topology = SubElement(self.arGrid,"Topology",
                                      {"Type":Xdmf_ElementTopology,
                                       "NumberOfElements":str(self.nElements_owned)})
                elements = SubElement(topology,"DataItem",
                                      {"Format":ar.dataItemFormat,
                                       "DataType":"Int",
                                       "Dimensions":"%i %i" % (self.nElements_owned,self.nNodes_element)})
                geometry = SubElement(self.arGrid,"Geometry",{"Type":"XYZ"})
                nodes    = SubElement(geometry,"DataItem",
                                      {"Format":ar.dataItemFormat,
                                       "DataType":"Float",
                                       "Precision":"8",
                                       "Dimensions":"%i %i" % (self.nNodes_global,3)})
                if ar.hdfFile is not None:
                    if ar.has_h5py:
                        elements.text = ar.hdfFilename+":/elements"+`ar.comm.rank()`+name+`tCount`
                        nodes.text = ar.hdfFilename+":/nodes"+`ar.comm.rank()`+name+`tCount`
                        if init or meshChanged:
                            ar.create_dataset_async('elements'+`ar.comm.rank()`+name+`tCount`,data=self.elementNodesArray[:self.nElements_owned])
                            ar.create_dataset_async('nodes'+`ar.comm.rank()`+name+`tCount`,data=self.nodeArray)
                    else:
                        elements.text = ar.hdfFilename+":/elements"+name+`tCount`
                        nodes.text = ar.hdfFilename+":/nodes"+name+`tCount`
                        if init or meshChanged:
                            ar.hdfFile.createArray("/",'elements'+name+`tCount`,self.elementNodesArray[:self.nElements_owned])
                            ar.hdfFile.createArray("/",'nodes'+name+`tCount`,self.nodeArray)
                else:
                    SubElement(elements,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/elements"+name+".txt"})
                    SubElement(nodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodes"+name+".txt"})
                    if init or meshChanged:
                        numpy.savetxt(ar.textDataDir+"/elements"+name+".txt",self.elementNodesArray[:self.nElements_owned],fmt='%d')
                        numpy.savetxt(ar.textDataDir+"/nodes"+name+".txt",self.nodeArray)
            #
            #element boundary topology and geometry
            #
            if EB:
                self.arEBGrid = SubElement(self.arEBGridCollection,"Grid",{"GridType":"Uniform"})
                self.arEBTime = SubElement(self.arEBGrid,"Time",{"Value":str(t),"Name":str(tCount)})
                Xdmf_ElementEBTopology = "Triangle" #cek hack
                ebtopology = SubElement(self.arEBGrid,"Topology",
                                    {"Type":Xdmf_ElementEBTopology,
                                     "NumberOfElements":str(self.nElementBoundaries_global)})
                ebelements = SubElement(ebtopology,"DataItem",
                                    {"Format":ar.dataItemFormat,
                                     "DataType":"Int",
                                     "Dimensions":"%i %i" % (self.nElementBoundaries_global,self.nNodes_elementBoundary)})
                ebgeometry = SubElement(self.arEBGrid,"Geometry",{"Type":"XYZ"})
                ebnodes    = SubElement(ebgeometry,"DataItem",
                                    {"Format":ar.dataItemFormat,
                                     "DataType":"Float",
                                     "Precision":"8",
                                     "Dimensions":"%i %i" % (self.nNodes_global,3)})
                if ar.hdfFile is not None:
                    if ar.has_h5py:
                        ebelements.text = ar.hdfFilename+":/elementBoundaries"+`ar.comm.rank()`+name+`tCount`
                        ebnodes.text = ar.hdfFilename+":/nodes"+`ar.comm.rank()`+name+`tCount`
                        if init or meshChanged:
                            ar.create_dataset_async('elementBoundaries'+`ar.comm.rank()`+name+`tCount`, data = self.elementBoundaryNodesArray)
                            #ar.create_dataset_async('nodes'+`ar.comm.rank()`+name+`tCount`, data = self.nodeArray)
                    else:
                        ebelements.text = ar.hdfFilename+":/elementBoundaries"+name+`tCount`
                        ebnodes.text = ar.hdfFilename+":/nodes"+name+`tCount`
                        if init or meshChanged:
                            ar.hdfFile.createArray("/",'elementBoundaries'+name+`tCount`,self.elementBoundaryNodesArray)
                            #ar.hdfFile.createArray("/",'nodes'+name+`tCount`,self.nodeArray)
                else:
                    SubElement(ebelements,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/elementBoundaries"+name+".txt"})
                    SubElement(ebnodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodes"+name+".txt"})
                    if init or meshChanged:
                        np.savetxt(ar.textDataDir+"/elementBoundaries"+name+".txt",self.elementBoundaryNodesArray,fmt='%d')

            # Add the local->global index maps for collect.py and for
            # reverse mapping in hotstarts from a global XDMF file.
            if self.globalMesh is not None and not ar.global_sync:
                nodeMapAtt = SubElement(self.arGrid,"Attribute",
                                        {"Name":"NodeMapL2G",
                                         "AttributeType":"Scalar",
                                         "Center":"Node"})
                nodeMap = SubElement(nodeMapAtt,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Int",
                                      "Precision":"4",
                                      "Dimensions":str(self.nNodes_global)})
                elemMapAtt = SubElement(self.arGrid,"Attribute",
                                        {"Name":"CellMapL2G",
                                         "AttributeType":"Scalar",
                                         "Center":"Cell"})
                elemMap = SubElement(elemMapAtt,"DataItem",
                                     {"Format":ar.dataItemFormat,
                                      "DataType":"Int",
                                      "Precision":"4",
                                      "Dimensions":str(self.nElements_owned)})

                if ar.hdfFile is not None:
                    if ar.has_h5py:
                        nodeMap.text = ar.hdfFilename+":/nodeMapL2G"+`ar.comm.rank()`+name+`tCount`
                        elemMap.text = ar.hdfFilename+":/cellMapL2G"+`ar.comm.rank()`+name+`tCount`
                        if init or meshChanged:
                            ar.create_dataset_async('nodeMapL2G'+`ar.comm.rank()`+name+`tCount`, data=self.globalMesh.nodeNumbering_subdomain2global)
                            ar.create_dataset_async('cellMapL2G'+`ar.comm.rank()`+name+`tCount`, data=self.globalMesh.elementNumbering_subdomain2global[:self.nElements_owned])
                    else:
                        nodeMap.text = ar.hdfFilename+":/nodeMapL2G"+name+`tCount`
                        elemMap.text = ar.hdfFilename+":/cellMapL2G"+name+`tCount`
                        if init or meshChanged:
                            ar.hdfFile.createArray("/",'nodeMapL2G'+name+`tCount`,self.globalMesh.nodeNumbering_subdomain2global)
                            ar.hdfFile.createArray("/",'cellMapL2G'+name+`tCount`,self.globalMesh.elementNumbering_subdomain2global[:self.nElements_owned])
                else:
                    SubElement(nodeMap,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodeMapL2G"+name+".txt"})
                    SubElement(nodeMap,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/cellMapL2G"+name+".txt"})
                    if init or meshChanged:
                        np.savetxt(ar.textDataDir+"/nodeMapL2G"+name+".txt",self.globalMesh.nodeNumbering_subdomain2global)
                        np.savetxt(ar.textDataDir+"/cellMapL2G"+name+".txt",self.globalMesh.elementNumbering_subdomain2global[:self.nElements_owned])
            #
            #material types
            #
            if ar.global_sync:
                nodeMaterialTypes = SubElement(self.arGrid,"Attribute",{"Name":"nodeMaterialTypes",
                                                                        "AttributeType":"Scalar",
                                                                        "Center":"Node"})
                nodeMaterialTypesValues = SubElement(nodeMaterialTypes,"DataItem",
                                                     {"Format":ar.dataItemFormat,
                                                      "DataType":"Int",
                                                      "Dimensions":"%i" % (self.globalMesh.nNodes_global,)})
                elementMaterialTypes = SubElement(self.arGrid,"Attribute",{"Name":"elementMaterialTypes",
                                                                           "AttributeType":"Scalar",
                                                                           "Center":"Cell"})
                elementMaterialTypesValues = SubElement(elementMaterialTypes,"DataItem",
                                                        {"Format":ar.dataItemFormat,
                                                         "DataType":"Int",
                                                         "Dimensions":"%i" % (self.globalMesh.nElements_global,)})
                if EB:
                    ebnodeMaterialTypes = SubElement(self.arEBGrid,"Attribute",{"Name":"ebnodeMaterialTypes",
                                                                          "AttributeType":"Scalar",
                                                                          "Center":"Node"})
                    ebnodeMaterialTypesValues = SubElement(ebnodeMaterialTypes,"DataItem",
                                                       {"Format":ar.dataItemFormat,
                                                        "DataType":"Int",
                                                        "Dimensions":"%i" % (self.globalMesh.nNodes_global,)})
                    elementBoundaryMaterialTypes = SubElement(self.arEBGrid,"Attribute",{"Name":"elementBoundaryMaterialTypes",
                                                                                             "AttributeType":"Scalar",
                                                                                             "Center":"Cell"})
                    elementBoundaryMaterialTypesValues = SubElement(elementBoundaryMaterialTypes,"DataItem",
                                                          {"Format":ar.dataItemFormat,
                                                           "DataType":"Int",
                                                           "Dimensions":"%i" % (self.globalMesh.nElementBoundaries_global,)})
                if ar.hdfFile is not None:
                    if ar.has_h5py:
                        nodeMaterialTypesValues.text = ar.hdfFilename+":/"+"nodeMaterialTypes"+"_t"+str(tCount)
                        ar.create_dataset_sync("nodeMaterialTypes"+"_t"+str(tCount), offsets=self.globalMesh.nodeOffsets_subdomain_owned, data=self.nodeMaterialTypes[:self.nNodes_owned])
                        elementMaterialTypesValues.text = ar.hdfFilename+":/"+"elementMaterialTypes"+"_t"+str(tCount)
                        ar.create_dataset_sync("elementMaterialTypes"+"_t"+str(tCount), offsets=self.globalMesh.elementOffsets_subdomain_owned, data=self.elementMaterialTypes[:self.nElements_owned])
                        if EB:
                            ebnodeMaterialTypesValues.text = ar.hdfFilename+":/"+"nodeMaterialTypes"+"_t"+str(tCount)
                            elementBoundaryMaterialTypesValues.text = ar.hdfFilename+":/"+"elementBoundaryMaterialTypes"+"_t"+str(tCount)
                            ar.create_dataset_sync("elementBoundaryMaterialTypes"+"_t"+str(tCount), offsets = self.globalMesh.elementBoundaryOffsets_subdomain_owned, data=self.elementBoundaryMaterialTypes[:self.nElementBoundaries_owned])
                    else:
                        assert False, "global_sync not supported  with pytables"
                else:
                    assert False, "global_sync  not  supported with text heavy data"
            else:
                nodeMaterialTypes = SubElement(self.arGrid,"Attribute",{"Name":"nodeMaterialTypes",
                                                                        "AttributeType":"Scalar",
                                                                        "Center":"Node"})
                nodeMaterialTypesValues = SubElement(nodeMaterialTypes,"DataItem",
                                                     {"Format":ar.dataItemFormat,
                                                      "DataType":"Int",
                                                      "Dimensions":"%i" % (self.nNodes_global,)})
                elementMaterialTypes = SubElement(self.arGrid,"Attribute",{"Name":"elementMaterialTypes",
                                                                           "AttributeType":"Scalar",
                                                                           "Center":"Cell"})
                elementMaterialTypesValues = SubElement(elementMaterialTypes,"DataItem",
                                                        {"Format":ar.dataItemFormat,
                                                         "DataType":"Int",
                                                         "Dimensions":"%i" % (self.nElements_owned,)})
                if EB:
                    ebnodeMaterialTypes = SubElement(self.arEBGrid,"Attribute",{"Name":"ebnodeMaterialTypes",
                                                                          "AttributeType":"Scalar",
                                                                          "Center":"Node"})
                    ebnodeMaterialTypesValues = SubElement(ebnodeMaterialTypes,"DataItem",
                                                       {"Format":ar.dataItemFormat,
                                                        "DataType":"Int",
                                                        "Dimensions":"%i" % (self.nNodes_global,)})
                    elementBoundaryMaterialTypes = SubElement(self.arEBGrid,"Attribute",{"Name":"elementBoundaryMaterialTypes",
                                                                                             "AttributeType":"Scalar",
                                                                                             "Center":"Cell"})
                    elementBoundaryMaterialTypesValues = SubElement(elementBoundaryMaterialTypes,"DataItem",
                                                          {"Format":ar.dataItemFormat,
                                                           "DataType":"Int",
                                                           "Dimensions":"%i" % (self.nElementBoundaries_global,)})
                if ar.hdfFile is not None:
                    if ar.has_h5py:
                        nodeMaterialTypesValues.text = ar.hdfFilename+":/"+"nodeMaterialTypes"+"_p"+`ar.comm.rank()`+"_t"+str(tCount)
                        ar.create_dataset_async("nodeMaterialTypes"+"_p"+`ar.comm.rank()`+"_t"+str(tCount), data=self.nodeMaterialTypes)
                        elementMaterialTypesValues.text = ar.hdfFilename+":/"+"elementMaterialTypes"+"_p"+`ar.comm.rank()`+"_t"+str(tCount)
                        ar.create_dataset_async("elementMaterialTypes"+"_p"+`ar.comm.rank()`+"_t"+str(tCount), data=self.elementMaterialTypes[:self.nElements_owned])
                        if EB:
                            ebnodeMaterialTypesValues.text = ar.hdfFilename+":/"+"nodeMaterialTypes"+"_p"+`ar.comm.rank()`+"_t"+str(tCount)
                            elementBoundaryMaterialTypesValues.text = ar.hdfFilename+":/"+"elementBoundaryMaterialTypes"+"_p"+`ar.comm.rank()`+"_t"+str(tCount)
                            ar.create_dataset_async("elementBoundaryMaterialTypes"+"_p"+`ar.comm.rank()`+"_t"+str(tCount), data=self.elementBoundaryMaterialTypes)
                    else:
                        nodeMaterialTypesValues.text = ar.hdfFilename+":/"+"nodeMaterialTypes"+str(tCount)
                        ar.hdfFile.createArray("/","nodeMaterialTypes"+str(tCount),self.nodeMaterialTypes)
                        elementMaterialTypesValues.text = ar.hdfFilename+":/"+"elementMaterialTypes"+str(tCount)
                        ar.hdfFile.createArray("/","elementMaterialTypes"+str(tCount),self.elementMaterialTypes[:self.nElements_owned])
                        if EB:
                            ebnodeMaterialTypesValues.text = ar.hdfFilename+":/"+"nodeMaterialTypes"+str(tCount)
                            #ar.hdfFile.createArray("/","nodeMaterialTypes"+str(tCount),self.nodeMaterialTypes)
                            elementBoundaryMaterialTypesValues.text = ar.hdfFilename+":/"+"elementBoundaryMaterialTypes"+str(tCount)
                            ar.hdfFile.createArray("/","elementBoundaryMaterialTypes"+str(tCount),self.elementBoundaryMaterialTypes)
                else:
                    numpy.savetxt(ar.textDataDir+"/"+"nodeMaterialTypes"+str(tCount)+".txt",self.nodeMaterialTypes)
                    SubElement(nodeMaterialTypesValues,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+"nodeMaterialTypes"+str(tCount)+".txt"})
                    numpy.savetxt(ar.textDataDir+"/"+"elementMaterialTypes"+str(tCount)+".txt",self.elementMaterialTypes[:self.nElements_owned])
                    SubElement(elementMaterialTypesValues,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+"elementMaterialTypes"+str(tCount)+".txt"})
            #done with material types
    def buildFromC(self,cmesh):
        import cmeshTools
        #
        logEvent(memory("buildFromC","MeshTools"),level=4)
        self.cmesh = cmesh
        (self.nElements_global,
         self.nNodes_global,
         self.nNodes_element,
         self.nNodes_elementBoundary,
         self.nElementBoundaries_element,
         self.nElementBoundaries_global,
         self.nInteriorElementBoundaries_global,
         self.nExteriorElementBoundaries_global,
         self.max_nElements_node,
         self.nEdges_global,
         self.max_nNodeNeighbors_node,
         self.elementNodesArray,
         self.nodeElementsArray,
         self.nodeElementOffsets,
         self.elementNeighborsArray,
         self.elementBoundariesArray,
         self.elementBoundaryNodesArray,
         self.elementBoundaryElementsArray,
         self.elementBoundaryLocalElementBoundariesArray,
         self.interiorElementBoundariesArray,
         self.exteriorElementBoundariesArray,
         self.edgeNodesArray,
         self.nodeStarArray,
         self.nodeStarOffsets,
         self.elementMaterialTypes,
         self.elementBoundaryMaterialTypes,
         self.nodeMaterialTypes,
         self.nodeArray,
         self.nx,self.ny, self.nz,      #NURBS
         self.px,self.py, self.pz,      #NURBS
         self.elementIJK, #NURBS
         self.weights,    #NURBS
         self.U_KNOT,     #NURBS
         self.V_KNOT,     #NURBS
         self.W_KNOT,     #NURBS
         self.elementDiametersArray,
         self.elementInnerDiametersArray,
         self.elementBoundaryDiametersArray,
         self.elementBarycentersArray,
         self.elementBoundaryBarycentersArray,
         self.nodeDiametersArray,
         self.nodeSupportArray,
         self.h,
         self.hMin,
         self.sigmaMax,
         self.volume) = cmeshTools.buildPythonMeshInterface(self.cmesh)
        self.hasGeometricInfo = True
        #default to single processor
        self.nNodes_owned = self.nNodes_global
        self.nElements_owned = self.nElements_global
        self.nElementBoundaries_owned = self.nElementBoundaries_global
        self.nEdges_owned = self.nEdges_global
        logEvent(memory("buildFromC","MeshTools"),level=4)
    def buildFromCNoArrays(self,cmesh):
        import cmeshTools
        #
        logEvent(memory("buildFromC","MeshTools"),level=4)
        self.cmesh = cmesh
        (self.nElements_global,
         self.nNodes_global,
         self.nNodes_element,
         self.nNodes_elementBoundary,
         self.nElementBoundaries_element,
         self.nElementBoundaries_global,
         self.nInteriorElementBoundaries_global,
         self.nExteriorElementBoundaries_global,
         self.max_nElements_node,
         self.nEdges_global,
         self.max_nNodeNeighbors_node,
         self.h,
         self.hMin,
         self.sigmaMax,
         self.volume) = cmeshTools.buildPythonMeshInterfaceNoArrays(self.cmesh)
        self.hasGeometricInfo = False
        logEvent(memory("buildFromCNoArrays","MeshTools"),level=4)
    def buildNodeStarArrays(self):
        if self.nodeStarArray is None:
            #cek old
            self.nodeStarList=[]
            for n in range(self.nNodes_global):
                self.nodeStarList.append([])
            for eNodes in self.edgeNodesArray:
                self.nodeStarList[eNodes[0]].append(eNodes[1])
                self.nodeStarList[eNodes[1]].append(eNodes[0])
            #cek new
            self.nodeStarOffsets = np.zeros((self.nNodes_global+1,),'i')
            lenNodeStarArray=0
            for nN in range(1,self.nNodes_global+1):
                self.nodeStarOffsets[nN] = self.nodeStarOffsets[nN-1] + len(self.nodeStarList[nN])
            self.nodeStarArray = np.array((self.nodeStarOffsets[-1],),'i')
            for nN in range(self.nNodes_global):
                for nN_star,offset in enumerate(range(self.nodeStarOffsets[nN],self.nodeStarOffsets[nN+1])):
                    self.nodeStarArray[offset] = self.nodeStarList[nN][nN_star]
            del self.nodeStarList
    def buildArraysFromLists(self):
        #nodes
        self.nNodes_global = len(self.nodeList)
        self.nodeArray = np.zeros((self.nNodes_global,3),'d')
        nodeElementsList=[]
        for nN,n in enumerate(self.nodeList):
            self.nodeArray[nN][:] = n.p
            nodeElementsList.append([])
        #elements
        self.nNodes_element = len(self.elementList[0].nodes)
        self.nElements_global = len(self.elementList)
        self.elementNodesArray = np.zeros((self.nElements_global,
                                                self.nNodes_element),
                                               'i')
        for en,e in enumerate(self.elementList):
            for nN_element,n in enumerate(e.nodes):
                self.elementNodesArray[en,nN_element]=n.N
                nodeElementsList[n.N].append(en)
        #elements per  node
        nodeElementsDict={}
        for eN in range(self.nElements_global):
            for nN_element in range(self.nNodes_element):
                nN = self.elementNodesArray[eN,nN_element]
                if nodeElementsDict.has_key(nN):
                    nodeElementsDict[nN].append(eN)
                else:
                    nodeElementsDict[nN] = [eN]
        self.max_nElements_node = max(len(nodeElementsDict[nN]) for  nN in range(self.nNodes_global))
        self.nElements_node = np.zeros((self.nNodes_global),'i')
        #mwf make a 1d array now
        #self.nodeElementsArray = np.zeros((self.nNodes_global,self.max_nElements_node),'i')
        self.nodeElementOffsets = np.zeros((self.nNodes_global+1,),'i')
        for nN,elementList in nodeElementsDict.iteritems():
            self.nElements_node[nN] = len(elementList)
            self.nodeElementOffsets[nN+1] = self.nodeElementOffsets[nN]+self.nElements_node[nN]
            #for eN_element,eN in enumerate(elementList):
            #    self.nodeElementsArray[nN,eN_element]=eN
        self.nodeElementsArray = np.zeros((self.nodeElementOffsets[self.nNodes_global],),'i')
        for nN,elementList in nodeElementsDict.iteritems():
            for eN_element,eN in enumerate(elementList):
                self.nodeElementsArray[self.nodeElementOffsets[nN]+eN_element]=eN
            #
        #
        #elementBoundariesArray
        self.nElementBoundaries_element = len(
            self.elementList[0].elementBoundaries)
        self.elementBoundariesArray = np.zeros(
            (self.nElements_global,self.nElementBoundaries_element),
            'i')
        #collect set of element boundaries while we're looping
        elementBoundaryNumbers=set()
        for eN,e in enumerate(self.elementList):
            for ebN_element,eb in enumerate(e.elementBoundaries):
                self.elementBoundariesArray[eN,ebN_element]=eb.N
                elementBoundaryNumbers.add(eb.N)
        self.nElementBoundaries_global=len(elementBoundaryNumbers)
        #elementBoundaryElementsArray
        self.elementBoundaryElementsArray=np.ones(
            (self.nElementBoundaries_global,2),'i')
        self.elementBoundaryElementsArray*=-1
        self.elementBoundaryLocalElementBoundariesArray=np.zeros(
            (self.nElementBoundaries_global,2),'i')
        elementBoundaryElementsCardArray =np.zeros(
            (self.nElementBoundaries_global),'i')
        for eN in range(self.nElements_global):
            for ebN_element in range(self.nElementBoundaries_element):
                ebN = self.elementBoundariesArray[eN,ebN_element]
                elementBoundaryElementsCardArray[ebN]+=1
                eN_boundaryElement=elementBoundaryElementsCardArray[ebN]-1
                self.elementBoundaryElementsArray[ebN,eN_boundaryElement]=eN
                self.elementBoundaryLocalElementBoundariesArray[ebN,eN_boundaryElement]=ebN_element
                if elementBoundaryElementsCardArray[ebN] > 2:
                    logEvent("WARNING, element neighbors of boundary element > 2")
                    elementBoundaryElementsCardArray[ebN]=2
        #interior and exterior
        self.nExteriorElementBoundaries_global=2*self.nElementBoundaries_global\
                                               - np.sum(
            elementBoundaryElementsCardArray)
        self.nInteriorElementBoundaries_global= self.nElementBoundaries_global-\
                                               self.nExteriorElementBoundaries_global
        self.exteriorElementBoundariesArray=np.zeros(
            (self.nExteriorElementBoundaries_global,),'i')
        self.interiorElementBoundariesArray=np.zeros(
            (self.nInteriorElementBoundaries_global,),'i')
        interior=0
        exterior=0
        for ebN in range(self.nElementBoundaries_global):
            if elementBoundaryElementsCardArray[ebN]==1:
                self.exteriorElementBoundariesArray[exterior]=ebN
                exterior+=1
            else:
                self.interiorElementBoundariesArray[interior]=ebN
                interior+=1
        del elementBoundaryElementsCardArray
        self.nNodes_elementBoundary = len(self.elementBoundaryList[0].nodes)
        self.elementBoundaryNodesArray = np.zeros((self.nElementBoundaries_global,
                                                        self.nNodes_elementBoundary),
                                                       'i')
        for ebN,eb in enumerate(self.elementBoundaryList):
            for nN_element,n in enumerate(eb.nodes):
                self.elementBoundaryNodesArray[ebN,nN_element]=n.N
        #element  neighbors
        self.elementNeighborsArray = np.zeros((self.nElements_global,self.nElementBoundaries_element),'i')
        for eN in range(self.nElements_global):
            for ebN_element in range(self.nElementBoundaries_element):
                ebN = self.elementBoundariesArray[eN,ebN_element]
                eN_left = self.elementBoundaryElementsArray[ebN,0]
                eN_right = self.elementBoundaryElementsArray[ebN,1]
                if eN == eN_left:
                    self.elementNeighborsArray[eN,ebN_element] = eN_right
                elif eN == eN_right:
                    self.elementNeighborsArray[eN,ebN_element] = eN_left
                else:
                    self.elementNeighborsArray[eN,ebN_element] = -1
        #edges
        self.edgeNodesArray = np.zeros(
            (len(self.edgeList),2),'i')
        for en,e in enumerate(self.edgeList):
            self.edgeNodesArray[en,0]=e.nodes[0].N
            self.edgeNodesArray[en,1]=e.nodes[1].N
        #geometric info
        self.computeGeometricInfo()
        self.elementDiametersArray = np.zeros((self.nElements_global,),'d')
        self.elementInnerDiametersArray = np.zeros((self.nElements_global,),'d')
        for en in range(self.nElements_global):
            self.elementDiametersArray[en] = self.elementList[en].diameter
            self.elementInnerDiametersArray[en]=self.elementList[en].innerDiameter
        self.elementBoundaryDiametersArray = np.zeros((self.nElementBoundaries_global,),'d')
        for eN,e in enumerate(self.elementList):
            for ebN_element,eb in enumerate(e.elementBoundaries):
                self.elementBoundaryDiametersArray[self.elementBoundariesArray[eN,ebN_element]] = eb.diameter
        self.elementMaterialTypes = np.zeros((self.nElements_global,),'i')
        self.elementBoundaryMaterialTypes = np.zeros((self.nElementBoundaries_global,),'i')
        self.nodeMaterialTypes = np.zeros((self.nNodes_global,),'i')
        #
        self.elementBarycentersArray         = np.zeros((self.nElements_global,3),'d')
        self.elementBoundaryBarycentersArray = np.zeros((self.nElementBoundaries_global,3),'d')
        for eN in range(self.nElements_global):
            self.elementBarycentersArray[eN,:] = 0.0
            for ebN in range(self.nNodes_element):
                self.elementBarycentersArray[eN,:] += self.nodeArray[self.elementNodesArray[eN,ebN],:]
            self.elementBarycentersArray[eN,:] /= float(self.nNodes_element)
        for ebN in range(self.nElementBoundaries_global):
            self.elementBoundaryBarycentersArray[ebN,:] = 0.0
            for nN in range(self.nNodes_elementBoundary):
                self.elementBoundaryBarycentersArray[ebN,:] += self.nodeArray[self.elementBoundaryNodesArray[ebN,nN],:]
            self.elementBoundaryBarycentersArray[ebN,:] /= float(self.nNodes_elementBoundary)
        #
        #now get rid of lists
        del self.nodeList
        del self.elementList
        del self.elementBoundaryList
        del self.edgeList
        #self.partitionMesh()
    def computeGeometricInfo(self):
        self.elementList[0].computeGeometricInfo()
        self.h=self.elementList[0].diameter
        self.hMin=self.h
        for e in self.elementList[1:]:
            e.computeGeometricInfo()
            self.h = max(self.h,e.diameter)
            self.hMin=min(self.hMin,e.diameter)
            for eb in e.elementBoundaries:
                e.computeGeometricInfo()
        self.hasGeometricInfo=True

    def buildMatlabMeshDataStructures(self,meshFileBase='meshMatlab',writeToFile=True):
        """
        build array data structures for matlab finite element mesh
        representation and write to a file to view and play with in
        matlatb. The current matlab support is mostly for 2d, but this
        will return basic arrays for 1d and 3d too

        in matlab can then print mesh with

        pdemesh(p,e,t)

        if one has pdetoolbox
        where

          p is the vertex or point matrix
          e is the edge matrix, and
          t is the element matrix

        e will be the elementBoundary matrix in 1d and 3d, but perhaps
        should remain the edge array?

        points matrix is [nd x num vertices]
          format :
             row 1 = x coord,
             row 2 = y coord for nodes in mesh
             row 3 = z coord for nodes in mesh ...

        edge matrix is [2*nd+3 x num faces]
          format:
             row 1  = start vertex number
             ...
             row nd   = end vertex number
             row nd+1 = start value in edge parameterization, should be 0
             row nd+2 = next value in edge parameterization, should be 1 or 2
             row nd+nd= end value in edge parameterization, should be 2 or 1
             row 2*nd+1 = global face id, base 1
             row 2*nd+2 = subdomain on left? always 1 for now
             row 2*nd+3 = subdomain on right? always 0 for now

        element matrix is [nd+2 x num elements]
            row 1 = vertex 1 global number
            row 2 = vertex 2 global number
            ...
            row nd+1 = vertex 3 global number
            row 4 = triangle subdomain number

         where 1,2,3 is a local counter clockwise numbering of vertices in
           triangle

        """
        matlabBase = 1
        nd = self.nNodes_element-1
        p = np.zeros((nd,self.nNodes_global),'d')
        e = np.zeros((2*nd+3,self.nElementBoundaries_global),'d')
        t = np.zeros((nd+2,self.nElements_global),'d')

        #load p,e,t and write file
        if writeToFile:
            mfile = open(meshFileBase+'.m','w')
        else:
            mfile = open('/dev/null','w')
        #
        if writeToFile:
            mfile.write('p = [ ... \n')
        for nN in range(self.nNodes_global):
            for I in range(nd):
                p[I,nN]=self.nodeArray[nN,I]
                if writeToFile:
                    mfile.write('%g ' % p[I,nN])
            mfile.write('\n')
        if writeToFile:
            mfile.write(']; \n')
            mfile.write("p = p\';\n")  #need transpose for matlab

        if writeToFile:
            mfile.write('e = [ ... \n')
        for ebN in range(self.nElementBoundaries_global):
            eN_left = self.elementBoundaryElementsArray[ebN,0]
            eN_right= self.elementBoundaryElementsArray[ebN,1]#-1 --> exterior
            for nN in range(self.nNodes_elementBoundary):
                e[nN,ebN]=self.elementBoundaryNodesArray[ebN,nN] + matlabBase #global node number of start node base 1
            #assume for now existing parameterization ok
            for nN in range(self.nNodes_elementBoundary):
                e[self.nNodes_elementBoundary+nN,ebN]=nN #edge param. is 0 to 1
            e[2*self.nNodes_elementBoundary+1,ebN] = ebN+matlabBase
            e[2*self.nNodes_elementBoundary+1,ebN] = self.elementMaterialTypes[eN_left] #subdomain to left
            if eN_right >= 0:
                e[2*self.nNodes_elementBoundary+2,ebN]= self.elementMaterialTypes[eN_right] #subdomain to right
            else:
                e[2*self.nNodes_elementBoundary+2,ebN]= -1
            if writeToFile:
                for i in range(e.shape[0]):
                    mfile.write(' %g ' % e[i,ebN])
                mfile.write(' \n ')
        if writeToFile:
            mfile.write(']; \n')
            mfile.write("e = e\';\n")  #need transpose for matlab

        #write triangles last
        if writeToFile:
            mfile.write('t = [ ... \n')
        for eN in range(self.nElements_global):
            for nN in range(self.nNodes_element):
                t[nN,eN]=self.elementNodesArray[eN,nN]+matlabBase    #global node number for vertex nN
            t[self.nNodes_element,eN]=self.elementMaterialTypes[eN]                     #subdomain id
            if writeToFile:
                for i in range(t.shape[0]):
                    mfile.write('%g ' % t[i,eN])
                mfile.write('\n')
        if writeToFile:
            mfile.write(']; \n');
            mfile.write("t = t\';\n") #need transpose for matlab


        mfile.close()
        return p,e,t

    def writeEdgesMatlab(self,filename):
        """store coordinates in files formatted for Matlab"""
        xfile=filename+'_x.grf'
        yfile=filename+'_y.grf'
        zfile=filename+'_z.grf'
        print 'Storing edge information in %s, %s, and %s' % \
              (xfile,yfile,zfile)
        xOut = open(xfile,'w')
        yOut = open(yfile,'w')
        zOut = open(zfile,'w')
        for edge in self.edgeList:
            xOut.write('%14.8e ' % edge.nodes[0].p[X] )
            yOut.write('%14.8e ' % edge.nodes[0].p[Y] )
            zOut.write('%14.8e ' % edge.nodes[0].p[Z] )
        xOut.write('\n')
        yOut.write('\n')
        zOut.write('\n')
        for edge in self.edgeList:
            xOut.write('%14.8e ' % edge.nodes[1].p[X])
            yOut.write('%14.8e ' % edge.nodes[1].p[Y])
            zOut.write('%14.8e ' % edge.nodes[1].p[Z])
        xOut.write('\n')
        yOut.write('\n')
        zOut.write('\n')
        xOut.close()
        yOut.close()
        zOut.close()

    def viewTetrahedraMatlab(self,filename):
        """plot the edges"""
        cmdfile = filename +'.m'
        xfile=filename+'_x.grf'
        yfile=filename+'_y.grf'
        zfile=filename+'_z.grf'
        xedges=filename+'_x'
        yedges=filename+'_y'
        zedges=filename+'_z'
        #the following is for debugging: plot each tet seperately
        nT = len(self.edgeList)/6
        plotcommand = "-r \"load " + xfile + \
                      ", load " + yfile + \
                      ", load " + zfile
        plots=''
        for i in range(nT):
            plots = plots + \
                    ", figure(" +`i+1`+")" \
                    ", axis([0 1 0 1 0 1]), plot3("+xedges+\
                    "(:,"+`i`+"*6+1:("+`i`+"+1)*6),"+yedges+\
                    "(:,"+`i`+"*6+1:("+`i`+"+1)*6),"+zedges+\
                    "(:,"+`i`+"*6+1:("+`i`+"+1)*6),\'b-\') "
        plotcommand = plotcommand + plots +'\"'
        cmdOut = open(cmdfile,'w')
        cmdOut.write(plotcommand)
        cmdOut.close()
        import os
        print 'Calling matlab to view mesh'
        os.execlp('matlab',
                  'matlab',
                  '-nodesktop',
                  '-nosplash',
                  '-r',
                  filename)

    def viewMeshMatlab(self,filename):
        """plot the edges"""
        cmdfile = filename +'.m'
        xfile=filename+'_x.grf'
        yfile=filename+'_y.grf'
        zfile=filename+'_z.grf'
        xedges=filename+'_x'
        yedges=filename+'_y'
        zedges=filename+'_z'
        plotcommand = "load " + xfile + \
                      ", load " + yfile + \
                      ", load " + zfile + \
                      ", figure " + \
                      ", axis([0 1 0 1 0 1]), plot3("+xedges+\
                      ","+yedges+\
                      ","+zedges+\
                      ",\'b-\')"
        print plotcommand
        cmdOut = open(cmdfile,'w')
        cmdOut.write(plotcommand)
        cmdOut.close()
        import os
        print 'Calling matlab to view mesh'
        os.execlp('matlab',
                  'matlab',
                  '-nodesktop',
                  '-nosplash',
                  '-r',
                  filename)
#          from os import popen
#          matlab = popen('matlab','w')
#          matlab.write(plotcommand+'\n')
#          matlab.flush()
#          raw_input('Please press return to continue...\n')

    def writeEdgesGnuplot(self,filename):
        """store coordinates in files formatted for Matlab"""
        datfile=filename+'.dat'
        print 'Storing edge information in %s' % datfile
        edgesOut = open(datfile,'w')
        for edge in self.edgeList:
            dataline = '%14.8e %14.8e %14.8e \n' % \
                       (edge.nodes[0].p[X],
                        edge.nodes[0].p[Y],
                        edge.nodes[0].p[Z])
            edgesOut.write(dataline)
            dataline = '%14.8e %14.8e %14.8e \n \n \n' % \
                       (edge.nodes[1].p[X],
                        edge.nodes[1].p[Y],
                        edge.nodes[1].p[Z])
            edgesOut.write(dataline)
        edgesOut.close()
    def writeEdgesGnuplot2(self,filename):
        """store coordinates in files formatted for Matlab"""
        datfile=filename+'.dat'
        print 'Storing edge information in %s' % datfile
        edgesOut = open(datfile,'w')
        for n0,n1 in self.edgeNodesArray:
            dataline = '%14.8e %14.8e %14.8e \n' % \
                       (self.nodeArray[n0][0],
                        self.nodeArray[n0][1],
                        self.nodeArray[n0][2])
            edgesOut.write(dataline)
            dataline = '%14.8e %14.8e %14.8e \n \n \n' % \
                       (self.nodeArray[n1][0],
                        self.nodeArray[n1][1],
                        self.nodeArray[n1][2])
            edgesOut.write(dataline)
        edgesOut.close()
    def viewMeshGnuplot(self,filename):
        cmdfile = filename +'.cmd'
        datfile = filename +'.dat'
        cmd = "set pointsize 2.5 \n set term x11 \n splot \'"+datfile+"\' with linespoints pointsize 2.5 pt 2\n"+\
              "set xlabel \'x\' \n set ylabel \'y\' \n set zlabel \'z\' \n "
        cmdOut = open(cmdfile,'w')
        cmdOut.write(cmd)
        cmdOut.close()
        from os import execlp
        print 'Calling gnuplot to view mesh'
        execlp('gnuplot','gnuplot',cmdfile,'-')
    def viewMeshGnuplotPipe(self,filename):
        cmdfile = filename +'.cmd'
        datfile = filename +'.dat'
        cmd = "set pointsize 1.5 \n set term x11 \n splot \'"+datfile+"\' with linespoints pointsize 2.5 pt 2 \n"+\
              "set xlabel \'x\' \n set ylabel \'y\' \n set zlabel \'z\' \n "
        cmdOut = open(cmdfile,'w')
        cmdOut.write(cmd)
        cmdOut.close()
        from os import execlp
        print 'Calling gnuplot to view mesh'
        from os import popen
        gnuplot = popen('gnuplot','w')
        gnuplot.write(cmd+'\n')
        gnuplot.flush()
        raw_input('Please press return to continue... \n')
    def viewMeshGnuplotPipePar(self,filenames):
        from os import popen
        gnuplot = popen('gnuplot','w')
        for i,filename in enumerate(filenames):
            cmdfile = filename +'.cmd'
            datfile = filename +'.dat'
            cmd = ("set term x11 %i \n splot \'" % (i,))+datfile+"\' with linespoints \n"+\
                  "set xlabel \'x\' \n set ylabel \'y\' \n set zlabel \'z\'"
            cmdOut = open(cmdfile,'w')
            cmdOut.write(cmd)
            cmdOut.close()
            from os import execlp
            print 'Calling gnuplot to view mesh'
            gnuplot.write(cmd+'\n')
            gnuplot.flush()
        raw_input('Please press return to continue... \n')

    def convertFromPUMI(self, PUMIMesh, faceList,regList, parallel=False, dim=3):
        import cmeshTools
        import MeshAdaptPUMI
        import flcbdfWrappers
        import Comm
        comm = Comm.get()
        self.cmesh = cmeshTools.CMesh()
        if parallel:
          self.subdomainMesh=self.__class__()
          self.subdomainMesh.globalMesh = self
          self.subdomainMesh.cmesh = cmeshTools.CMesh()
          PUMIMesh.constructFromParallelPUMIMesh(self.cmesh,
              self.subdomainMesh.cmesh)
          if(PUMIMesh.isReconstructed()):
            logEvent("Material arrays updating based on reconstructed model.\n")
            PUMIMesh.updateMaterialArrays(self.subdomainMesh.cmesh);
          else:
              logEvent("Material arrays updating based on geometric model.\n")
              for i in range(len(faceList)):
                for j in range(len(faceList[i])):
                  PUMIMesh.updateMaterialArrays(self.subdomainMesh.cmesh,(dim-1), i+1,
                      faceList[i][j])
              for i in range(len(regList)):
                for j in range(len(regList[i])):
                  PUMIMesh.updateMaterialArrays(self.subdomainMesh.cmesh,dim, i+1, regList[i][j])
          if dim == 3:
            cmeshTools.allocateGeometricInfo_tetrahedron(self.subdomainMesh.cmesh)
            cmeshTools.computeGeometricInfo_tetrahedron(self.subdomainMesh.cmesh)
          if dim == 2:
            cmeshTools.allocateGeometricInfo_triangle(self.subdomainMesh.cmesh)
            cmeshTools.computeGeometricInfo_triangle(self.subdomainMesh.cmesh)
          self.buildFromCNoArrays(self.cmesh)
          (self.elementOffsets_subdomain_owned,
           self.elementNumbering_subdomain2global,
           self.nodeOffsets_subdomain_owned,
           self.nodeNumbering_subdomain2global,
           self.elementBoundaryOffsets_subdomain_owned,
           self.elementBoundaryNumbering_subdomain2global,
           self.edgeOffsets_subdomain_owned,
           self.edgeNumbering_subdomain2global) = (
              flcbdfWrappers.convertPUMIPartitionToPython(self.cmesh,
                  self.subdomainMesh.cmesh))
          self.subdomainMesh.buildFromC(self.subdomainMesh.cmesh)
          self.subdomainMesh.nElements_owned = (
              self.elementOffsets_subdomain_owned[comm.rank()+1] -
              self.elementOffsets_subdomain_owned[comm.rank()])
          self.subdomainMesh.nNodes_owned = (
              self.nodeOffsets_subdomain_owned[comm.rank()+1] -
              self.nodeOffsets_subdomain_owned[comm.rank()])
          self.subdomainMesh.nElementBoundaries_owned = (
              self.elementBoundaryOffsets_subdomain_owned[comm.rank()+1] -
              self.elementBoundaryOffsets_subdomain_owned[comm.rank()])
          self.subdomainMesh.nEdges_owned = (
              self.edgeOffsets_subdomain_owned[comm.rank()+1] -
              self.edgeOffsets_subdomain_owned[comm.rank()])
          comm.barrier()
          par_nodeDiametersArray = (
              ParVec_petsc4py(self.subdomainMesh.nodeDiametersArray,
                              bs=1,
                              n=self.subdomainMesh.nNodes_owned,
                              N=self.nNodes_global,
                              nghosts = self.subdomainMesh.nNodes_global -
                                        self.subdomainMesh.nNodes_owned,
                              subdomain2global = 
                                  self.nodeNumbering_subdomain2global))
          par_nodeDiametersArray.scatter_forward_insert()
          comm.barrier()
        else:
          PUMIMesh.constructFromSerialPUMIMesh(self.cmesh)
          if(PUMIMesh.isReconstructed()):
            PUMIMesh.updateMaterialArrays(self.cmesh);
          else:
              for i in range(len(faceList)):
                for j in range(len(faceList[i])):
                  PUMIMesh.updateMaterialArrays(self.cmesh,(dim-1), i+1, faceList[i][j])
              for i in range(len(regList)):
                for j in range(len(regList[i])):
                  PUMIMesh.updateMaterialArrays(self.cmesh,dim, i+1, regList[i][j])
          if dim == 3:
            cmeshTools.allocateGeometricInfo_tetrahedron(self.cmesh)
            cmeshTools.computeGeometricInfo_tetrahedron(self.cmesh)
          if dim == 2:
            cmeshTools.allocateGeometricInfo_triangle(self.cmesh)
            cmeshTools.computeGeometricInfo_triangle(self.cmesh)
          self.buildFromC(self.cmesh)
        logEvent("meshInfo says : \n"+self.meshInfo())

class MultilevelMesh(Mesh):
    """A hierchical multilevel mesh"""
    def __init__(self,levels=1):
        self.meshList=[]
        self.elementParents=None
    def buildFromC(self,cmultilevelMesh):
        import cmeshTools
        self.cmultilevelMesh = cmultilevelMesh
        (self.nLevels,
         self.cmeshList,
         self.elementParentsArrayList,
         self.elementChildrenArrayList,
         self.elementChildrenOffsetsList) = cmeshTools.buildPythonMultilevelMeshInterface(cmultilevelMesh)
    def refine(self):
        pass
    def locallyRefine(self,elementTagArray):
        pass
    def buildArrayLists(self):
        self.nLevels = len(self.meshList)
        self.calculateElementParents()
        self.elementParentsArrayList=[[]]
        self.elementChildrenArrayList=[]
        self.elementChildrenOffsetsList=[]
        for l in range(1,self.nLevels):
            self.elementParentsArrayList.append(self.elementParents[l])
            len_children=0
            for children in self.elementChildren[l-1].values():
                len_children += len(children)
            self.elementChildrenArrayList.append(np.zeros((len_children,),'i'))
            self.elementChildrenOffsetsList.append(np.zeros((self.meshList[l-1].nElements_global+1,),'i'))
            index=0
            for eN_p,children in enumerate(self.elementChildren[l-1].values()):
                self.elementChildrenOffsetsList[l-1][eN_p] = index
                for ec in children:
                    self.elementChildrenArrayList[l-1][index] = ec.N
                    index += 1
            self.elementChildrenOffsetsList[l-1][-1] = index
    def calculateElementParents(self,recalculate=False):
        """
        get array elementParents[l,e] = e_c, where element e_c is the parent of element e
            elementParents[0,:] = -1
        """
        if (self.elementParents is None or recalculate):
            self.elementParents = {}
            nLevels = len(self.meshList)
            for l in range(nLevels):
                nE   = self.meshList[l].nElements_global
                self.elementParents[l] = np.ones((nE,),'i')
                self.elementParents[l][:]=-1
            for l in range(0,nLevels-1):
                nEc = self.meshList[l].nElements_global
                for ec in range(nEc):
                    for ef in self.elementChildren[l][ec]:
                        #print """l=%s ec= %s ef.N= %s """ % (l,ec,ef.N)
                        self.elementParents[l+1][ef.N] = ec
                    #ef
                #ec
            #l

class PointMesh(Mesh):
    #Elements=Nodes
    """
    0D mesh
    """
    def __init__(self,points):
        self.nodeArray=points
        self.nNodes_global = points.shape[0]
        self.elementNodesArray=np.arange(self.nNodes_global,dtype='i')
        self.nElements_global = self.nNodes_global

MX=0
MY=1
MZ=2
I=0
J=1
K=1

class EdgeGrid(Mesh):
    """A 1D regular grid on an interval"""
    def __init__(self,nx=2,Lx=1.0):
        Mesh.__init__(self)
        #dimensions and ranges
        self.nx=nx
        self.ex=nx-1
        self.nRange_x = range(nx)
        self.eRange_x = range(self.ex)
        #lengths
        self.Lx=Lx
        self.dx = Lx/self.ex
        #node coordinates
        self.nodeGridArray = np.zeros((self.nx,3),'d')
        for i in self.nRange_x:
            self.nodeGridArray[i,MX] = i*self.dx
        #edge node numbers
        self.edgeNodesArray=np.zeros((self.ex,2),'i')
        #try to do this like we'll do 2d and 3d
        #edge nodes
        en=2
        edgeNodeNumbers = np.zeros((en,),'i')
        #reference edge
        eI=1
        refEdge_nodeIndeces = [-eI,eI]
        refEdge_NodeDict={}
        for rn,rnii in enumerate(refEdge_nodeIndeces):
            refEdge_NodeDict[rnii] = rn
        for i in self.eRange_x:
            #edge number
            eN=i
            #fine grid index of edge
            ii = 2*i + 1
            #fine grid index of edge nodes
            for rn,rnii in enumerate(refEdge_nodeIndeces):
                nii = rnii + ii
                edgeNodeNumbers[rn]=nii/2
            self.edgeNodesArray[eN,:]=edgeNodeNumbers
        #Mesh interface
        self.nNodes_global=self.nx
        self.nEdges_global=self.ex
        self.nElements_global=self.ex
        self.nElementBoundaries_global=self.nx
        self.nodeArray=self.nodeGridArray
        self.elementNodesArray=self.edgeNodesArray
        self.elementBoundariesArray=self.nodeArray
        self.boundaryMesh=PointMesh(np.array([self.nodeArray[0],
                                              self.nodeArray[-1]],dtype='d'))

class QuadrilateralGrid(Mesh):
    """A 2D regular grid of quadrilateral cells"""
    def __init__(self,nx=2,ny=2,Lx=1.0,Ly=1.0):
        Mesh.__init__(self)
        #nodes
        self.nx=nx
        self.ny=ny
        self.nxy=nx*ny
        #edges
        self.eXx=nx-1
        self.eXy=ny
        self.eXxy=self.eXx*self.eXy
        self.eYx=nx
        self.eYy=ny-1
        self.eYxy = self.eYx*self.eYy
        self.eXYx = self.eXx + self.eYx
        self.eXYy = self.eXy + self.eYy
        self.eXYxy = self.eXxy + self.eYxy
        #quads
        self.qx = nx-1
        self.qy = ny-1
        self.qxy = self.qx*self.qy
        #ranges
        self.nRange_x = range(self.nx)
        self.nRange_y = range(self.ny)
        self.qRange_x = range(self.qx)
        self.qRange_y = range(self.qx)
        #lengths
        self.Lx=Lx
        self.Ly=Ly
        self.dx = Lx/self.eXx
        self.dy = Ly/self.eYy
        #node coordinates
        self.nodeGridArray=np.zeros((nx,ny,3),'d')
        for i in self.nRange_x:
            for j in self.nRange_y:
                self.nodeGridArray[i,j,MX]=i*self.dx
                self.nodeGridArray[i,j,MY]=j*self.dy
        #edge node numbers
        en=2
        edgeNodeNumbers = np.zeros((en,),'i')
        self.edgeNodesArray=np.zeros((self.eXYxy,en),'i')
        #quad node numbers
        qn=4
        quadNodeNumbers = np.zeros((qn,),'i')
        self.quadrilateralNodesArray=np.zeros((self.qxy,qn),'i')
        #quad edge numbers
        qe=4
        quadEdgeNumbers = np.zeros((qe,),'i')
        self.quadrilateralEdgesArray=np.zeros((self.qxy,qe),'i')
        #reference quad
        refQuad_NodeIndeces = [(-1,-1),
                               (-1, 1),
                               ( 1,-1),
                               ( 1, 1)]
        refQuad_NodeIndeces.sort()
        #a map between reference node indeces and numbers
        refQuad_NodeDict={}
        for rn,rniijj in enumerate(refQuad_NodeIndeces):
            refQuad_NodeDict[rniijj]=rn
        refQuad_EdgeIndeces = [(-1,0),
                               ( 0,-1),
                               ( 0, 1),
                               ( 1, 0)]
        refQuad_EdgeIndeces.sort()
        refQuad_EdgeNodes=[]
        #use the map between indeces and numbers to
        #map edge indeces to the edge's node numbers
        for reiijj in refQuad_EdgeIndeces:
            if reiijj[I] == 0:
                refQuad_EdgeNodes.append([
                    refQuad_NodeDict[(-1,reiijj[J])],
                    refQuad_NodeDict[( 1,reiijj[J])]])
            else:
                refQuad_EdgeNodes.append([
                    refQuad_NodeDict[(reiijj[I],-1)],
                    refQuad_NodeDict[(reiijj[I], 1)]])
        for i in self.qRange_x:
            for j in self.qRange_y:
                #quad number
                qN = i*self.qy + j
                #fine grid indeces of quad
                ii = 2*i + 1
                jj = 2*j + 1
                #nodes
                for rn,rniijj in enumerate(refQuad_NodeIndeces):
                    nii = rniijj[I] + ii
                    njj = rniijj[J] + jj
                    nN = (nii/2)*self.ny + njj/2
                    quadNodeNumbers[rn]=nN
                self.quadrilateralNodesArray[qN][:]=quadNodeNumbers
                #edges
                for re,reiijj in enumerate(refQuad_EdgeIndeces):
                    eii = reiijj[I] + ii
                    ejj = reiijj[J] + jj
                    eN = (eii/2)*self.eXYy + (eii%2)*self.eYy + ejj/2
                    quadEdgeNumbers[re]=eN
                    #nodes
                    for n,rn in enumerate(refQuad_EdgeNodes[re]):
                        self.edgeNodesArray[eN][n] = quadNodeNumbers[rn]
                self.quadrilateralEdgesArray[qN][:]=quadEdgeNumbers
        #Mesh interface (dimensions)
        self.nNodes_global=self.nxy
        self.nEdges_global=self.eXYxy
        self.nElements_global=self.qxy
        self.nElementBoundaries_global=self.eXYxy
        self.nodeArray=np.reshape(self.nodeGridArray,(self.nxy,3))
        self.elementNodesArray=self.quadrilateralNodesArray
        self.elementBoundariesArray=self.edgeNodesArray
        #todo extract boundary mesh

class RectangularGrid(Mesh):
    """A regular partition into rectangles.

    Nodes, edges, and faces can be indexed by (i,j,k) as follows.
    The edges and faces are divided according to orientation (i.e. x-edge...).
    An (i,j,k) index is associated with the type of edge or face
    having node (i,j,k) as the first node in a lexicographically sorted
    list of nodes corresponding to the edge or face."""
    def __init__(self,nx=1,ny=1,nz=1,Lx=1.0,Ly=1.0,Lz=1.0):
        Mesh.__init__(self)
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        #nodes
        self.nx=nx
        self.ny=ny
        self.nz=nz
        self.nxy = nx*ny;
        self.nxyz = nx*ny*nz;

        #edges
        self.nXex=nx-1 #number of x-edges in the x dimension
        self.nXey=ny
        self.nXez=nz

        self.nYex=nx
        self.nYey=ny-1
        self.nYez=nz

        self.nZex=nx
        self.nZey=ny
        self.nZez=nz-1

        #number of edges of all types associated with a row of nodes
        self.nXYZex = self.nXex + self.nYex + self.nZex

        #number of edges associated with an xy plane of nodes
        self.nXexy=self.nXex*self.nXey
        self.nYexy=self.nYex*self.nYey
        self.nZexy=self.nZex*self.nZey
        self.nXYZexy = self.nXexy + self.nYexy + self.nZexy

        #number of edges of each type in the grid
        self.nXexyz = self.nXexy*self.nXez
        self.nYexyz = self.nYexy*self.nYez
        self.nZexyz = self.nZexy*self.nZez

        #total number of edges
        self.nXYZexyz = self.nXexyz + self.nYexyz + self.nZexyz

        #quadrilaterals
        self.nXYhx=nx-1 #number of XY quadrilaterals in x-dimension
        self.nXYhy=ny-1
        self.nXYhz=nz

        self.nXZhx=nx-1
        self.nXZhy=ny
        self.nXZhz=nz-1

        self.nYZhx=nx
        self.nYZhy=ny-1
        self.nYZhz=nz-1

        #number of quadrilaterals of all types associate with a row of nodes
        self.nXY_XZ_YZhx =self.nXYhx + self.nXZhx + self.nYZhx

        #number of quadrilaterals associated with an xy plane of nodes
        self.nXYhxy=self.nXYhx*self.nXYhy
        self.nXZhxy=self.nXZhx*self.nXZhy
        self.nYZhxy=self.nYZhx*self.nYZhy
        self.nXY_XZ_YZhxy =self.nXYhxy + self.nXZhxy + self.nYZhxy

        #number of quadrilaterals of each type in the grid
        self.nXYhxyz = self.nXYhxy*self.nXYhz
        self.nXZhxyz = self.nXZhxy*self.nXZhz
        self.nYZhxyz = self.nYZhxy*self.nYZhz

        #total number of quadrilaterals
        self.nXY_XZ_YZhxyz =self.nXYhxyz + self.nXZhxyz + self.nYZhxyz

        #hexahedra
        self.nHx=nx-1
        self.nHy=ny-1
        self.nHz=nz-1
        self.nHxy = self.nHx*self.nHy
        self.nHxyz = self.nHxy*self.nHz

        #encode red and black
        self.black=0
        self.red=1

        #dimensions of hexahedra
        if self.nHx>0:
            hx = float(Lx)/(nx-1)
        else:
            hx = 1.0
        if self.nHy>0:
            hy = float(Ly)/(ny-1)
        else:
            hy=1.0
        if self.nHz>0:
            hz = float(Lz)/(nz-1)
        else:
            hz=1.0
        self.nodeDict={}
        self.xedgeDict={}
        self.yedgeDict={}
        self.zedgeDict={}
        self.xedgeList=[]
        self.yedgeList=[]
        self.zedgeList=[]
        self.XYQuadrilateralDict={}
        self.XZQuadrilateralDict={}
        self.YZQuadrilateralDict={}
        self.XYQuadrilateralList=[]
        self.XZQuadrilateralList=[]
        self.YZQuadrilateralList=[]
        self.hexahedronDict={}
        self.hexahedronList=[]
        self.nodeList=[]
        for k in range(self.nz):
            for j in range(self.ny):
                for i in range(self.nx):
                    n = self.getNodeNumber(i,j,k)
                    x = i*hx
                    y = j*hy
                    z = k*hz
                    self.nodeDict[(i,j,k)]=Node(n,x,y,z)
                    self.nodeList.append(self.nodeDict[(i,j,k)])
        for k in range(self.nXez):
            for j in range(self.nXey):
                for i in range(self.nXex):
                    en = self.getXEdgeNumber(i,j,k)
                    self.xedgeDict[(i,j,k)] = Edge(en,
                                                   [self.getNode(i,j,k),
                                                    self.getNode(i+1,j,k)])
                    self.xedgeList.append(self.xedgeDict[(i,j,k)])
        for k in range(self.nYez):
            for j in range(self.nYey):
                for i in range(self.nYex):
                    en = self.getYEdgeNumber(i,j,k)
                    self.yedgeDict[(i,j,k)] = Edge(en,
                                                   [self.getNode(i,j,k),
                                                    self.getNode(i,j+1,k)])
                    self.yedgeList.append(self.yedgeDict[(i,j,k)])
        for k in range(self.nZez):
            for j in range(self.nZey):
                for i in range(self.nZex):
                    en = self.getZEdgeNumber(i,j,k)
                    self.zedgeDict[(i,j,k)] = Edge(en,
                                                   [self.getNode(i,j,k),
                                                    self.getNode(i,j,k+1)])
                    self.zedgeList.append(self.zedgeDict[(i,j,k)])
        for k in range(self.nXYhz):
            for j in range(self.nXYhy):
                for i in range(self.nXYhx):
                    qn = self.getXYQuadrilateralNumber(i,j,k)
                    edges = [self.getXEdge(i,j,k),
                             self.getXEdge(i,j+1,k),
                             self.getYEdge(i,j,k),
                             self.getYEdge(i+1,j,k)]
                    self.XYQuadrilateralDict[(i,j,k)] = Quadrilateral(qn,edges)
                    self.XYQuadrilateralList.append(
                        self.XYQuadrilateralDict[(i,j,k)])
        for k in range(self.nXZhz):
            for j in range(self.nXZhy):
                for i in range(self.nXZhx):
                    qn = self.getXZQuadrilateralNumber(i,j,k)
                    edges = [self.getXEdge(i,j,k),
                             self.getXEdge(i,j,k+1),
                             self.getZEdge(i,j,k),
                             self.getZEdge(i+1,j,k)]
                    self.XZQuadrilateralDict[(i,j,k)] = Quadrilateral(qn,edges)
                    self.XZQuadrilateralList.append(
                        self.XZQuadrilateralDict[(i,j,k)])
        for k in range(self.nYZhz):
            for j in range(self.nYZhy):
                for i in range(self.nYZhx):
                    qn = self.getYZQuadrilateralNumber(i,j,k)
                    edges = [self.getYEdge(i,j,k),
                             self.getYEdge(i,j,k+1),
                             self.getZEdge(i,j,k),
                             self.getZEdge(i,j+1,k)]
                    self.YZQuadrilateralDict[(i,j,k)] = Quadrilateral(qn,edges)
                    self.YZQuadrilateralList.append(
                        self.YZQuadrilateralDict[(i,j,k)])
        for  k in range(self.nHz):
            for j in range(self.nHy):
                for i in range(self.nHx):
                    Hn = self.getHexahedronNumber(i,j,k)
                    quadrilaterals = [self.getXYQuadrilateral(i,j,k),
                                      self.getXYQuadrilateral(i,j,k+1),
                                      self.getXZQuadrilateral(i,j,k),
                                      self.getXZQuadrilateral(i,j+1,k),
                                      self.getYZQuadrilateral(i,j,k),
                                      self.getYZQuadrilateral(i+1,j,k)]
                    self.hexahedronDict[(i,j,k)] = Hexahedron(Hn,
                                                              quadrilaterals)
                    self.hexahedronList.append(self.hexahedronDict[(i,j,k)])
        #build lists for mesh base class
        self.edgeList = self.xedgeList + \
                        self.yedgeList + \
                        self.zedgeList
        #figure out if this is a 1D,2D, or 3D grid
        if self.nz > 1:
            self.elementList = self.hexahedronList
            self.elementDict = self.hexahedronDict
        elif self.ny > 1:
            self.elementList = self.XYQuadrilateralList
            self.elementDict = self.XYQuadrilateralDict
        else:
            self.elementList = self.xedgeList
            self.elementDict = self.xedgeDict
        #self.buildArraysFromLists()
        #todo: extract boundary mesh

    def getNodeNumber(self,i,j,k):
        return i + j*self.nx + k*self.nxy

    def getNode(self,i,j,k):
        return self.nodeDict[(i,j,k)]

    def getXEdgeNumber(self,ie,je,ke):
        return ie + je*self.nXex + ke*self.nXexy

    def getYEdgeNumber(self,ie,je,ke):
        return ie + je*self.nYex + ke*self.nYexy

    def getZEdgeNumber(self,ie,je,ke):
        return ie + je*self.nZex + ke*self.nZexy

    def getXEdge(self,ie,je,ke):
        return self.xedgeDict[(ie,je,ke)]

    def getYEdge(self,ie,je,ke):
        return self.yedgeDict[(ie,je,ke)]

    def getZEdge(self,ie,je,ke):
        return self.zedgeDict[(ie,je,ke)]

    def getXYQuadrilateralNumber(self,ih,jh,kh):
        return ih + jh*self.nXYhx + kh*self.nXYhxy

    def getXZQuadrilateralNumber(self,ih,jh,kh):
        return ih + jh*self.nXZhx + kh*self.nXZhxy

    def getYZQuadrilateralNumber(self,ih,jh,kh):
        return ih + jh*self.nYZhx + kh*self.nYZhxy

    def getXYQuadrilateral(self,ih,jh,kh):
        return self.XYQuadrilateralDict[(ih,jh,kh)]

    def getXZQuadrilateral(self,ih,jh,kh):
        return self.XZQuadrilateralDict[(ih,jh,kh)]

    def getYZQuadrilateral(self,ih,jh,kh):
        return self.YZQuadrilateralDict[(ih,jh,kh)]

    def getHexahedronNumber(self,iH,jH,kH):
        return iH + jH*self.nHx + kH*self.nHxy

    def getHexahedron(self,iH,jH,kH):
        return self.hexahedronDict[(iH,jH,kH)]

    def getColor(self,i,j,k):
        return (i%2 + j%2 + k%2)%2

    def refine(self,oldMesh,refineFactorX=2,refineFactorY=2,refineFactorZ=2):
        NX = oldMesh.nx
        NY = oldMesh.ny
        NZ = oldMesh.nz
        if NX > 1:
            NX = (NX-1)*refineFactorX + 1
        else:
            refineFactorX=1
        if NY > 1:
            NY = (NY-1)*refineFactorY + 1
        else:
            refineFactorY=1
        if NZ > 1:
            NZ = (NZ-1)*refineFactorZ + 1
        else:
            refineFactorZ=1
        RectangularGrid.__init__(self,NX,NY,NZ,
                                 oldMesh.Lx,oldMesh.Ly,oldMesh.Lz)
        childrenDict={}
        for IJK,e in oldMesh.elementDict.iteritems():
            I = IJK[0]
            J = IJK[1]
            K = IJK[2]
            childrenDict[e.N]=[]
            for xOffset in range(refineFactorX):
                for yOffset in range(refineFactorY):
                    for zOffset in range(refineFactorZ):
                        i = I*refineFactorX + xOffset
                        j = J*refineFactorY + yOffset
                        k = K*refineFactorZ + zOffset
                        childrenDict[e.N].append(self.elementDict[(i,j,k)])
        return childrenDict

class MultilevelRectangularGrid(MultilevelMesh):
    """A hierarchical multilevel grid"""
    def __init__(self,levels,nx,ny=1,nz=1,
                 Lx=1.0,Ly=1.0,Lz=1.0,
                 refinementLevels=1):
        MultilevelMesh.__init__(self)
        self.refineFactorList=[EVec(0,0,0)]
        self.meshList.append(RectangularGrid(nx,ny,nz,Lx,Ly,Lz))
        self.elementChildren = []
        logEvent(self.meshList[0].meshInfo())
        for l in range(1,refinementLevels+1):
            self.refine()
            logEvent(self.meshList[-1].meshInfo())

    def refine():
        self.meshList.append(RectangularMesh())
        childrenDict = self.meshList[-1].refine(self.meshList[-2])
        self.elementChildren.append(childrenDict)

class TetrahedralMesh(Mesh):
    """A mesh of tetrahedra.

    The nodes, edges, triangles, and tetrahedra are indexed by their
    node tuples. The corresponding lists are derived from the dictionaries, and
    sorted lexicographically. The global node numbers are redefined to
    give a lexicographic ordering.

    The mesh can be generated from a rectangular grid and refined using either
    4T or Freudenthal-Bey global refinement.

    Attributes
    ----------
    elementNodesArray : array_like
        A list of lists storing the node values associated with each element 
        in the triangulation.  The first index refers to the element number,
        while the second index refers to the global node value.
    nodeArray : array_like
        A list of lists storing node coordinates.  The first index referes
        to the global node number, while the second index refers to the x, y
        and z coordinates of the node respectively.
    """

    def __init__(self):
        Mesh.__init__(self)
        self.nodeDict={}
        self.edgeDict={}
        self.triangleDict={}
        self.triangleList=[]
        self.tetrahedronDict={}
        self.tetrahedronList=[]
        self.oldToNewNode=[]
        self.boundaryMesh=TriangularMesh()
    def meshType(self):
        return 'simplex'
    def computeGeometricInfo(self):
        import cmeshTools
        cmeshTools.computeGeometricInfo_tetrahedron(self.cmesh)
    def generateTetrahedralMeshFromRectangularGrid(self,nx,ny,nz,Lx,Ly,Lz):
        import cmeshTools
        self.cmesh = cmeshTools.CMesh()
        cmeshTools.generateTetrahedralMeshFromRectangularGrid(nx,ny,nz,Lx,Ly,Lz,self.cmesh)
        cmeshTools.allocateGeometricInfo_tetrahedron(self.cmesh)
        cmeshTools.computeGeometricInfo_tetrahedron(self.cmesh)
        self.buildFromC(self.cmesh)
        cmeshTools.writeTetgenFiles(self.cmesh,"tetgen",1)
    def rectangularToTetrahedral6T(self,grid):
        #copy the nodes from the rectangular mesh
        #I want to be able to renumber later without
        #changing the grid nodes, so I do deep copies here
        self.nodeList = [Node(n.N,n.p[X],n.p[Y],n.p[Z]) for n in grid.nodeList]
        self.nodeDict = dict([(n,n) for n in self.nodeList])
        for i in range(grid.nHx):
            for j in range(grid.nHy):
                for k in range(grid.nHz):
                    #associate the element (i,j,k) with the
                    #left, front, bottom node
                    #do a top down numbering to match Ong's dissertation
                    n1 = self.nodeList[grid.getNodeNumber(i,j,k+1)]
                    n2 = self.nodeList[grid.getNodeNumber(i,j+1,k+1)]
                    n3 = self.nodeList[grid.getNodeNumber(i+1,j+1,k+1)]
                    n4 = self.nodeList[grid.getNodeNumber(i+1,j,k+1)]
                    n5 = self.nodeList[grid.getNodeNumber(i,j,k)]
                    n6 = self.nodeList[grid.getNodeNumber(i,j+1,k)]
                    n7 = self.nodeList[grid.getNodeNumber(i+1,j+1,k)]
                    n8 = self.nodeList[grid.getNodeNumber(i+1,j,k)]
                    self.newTetrahedron(nodes=[n1,n2,n3,n6])
                    self.newTetrahedron(nodes=[n1,n3,n5,n6])
                    self.newTetrahedron(nodes=[n3,n5,n6,n7])
                    self.newTetrahedron(nodes=[n1,n3,n4,n5])
                    self.newTetrahedron(nodes=[n3,n4,n5,n7])
                    self.newTetrahedron(nodes=[n4,n5,n7,n8])
        self.finalize()

    def rectangularToTetrahedral5T(self,grid):
        #copy the nodes from the rectangular mesh
        #I want to be able to renumber later without
        #changing the grid nodes, so I do deep copies here
        self.nodeList = [Node(n.N,n.p[X],n.p[Y],n.p[Z]) for n in grid.nodeList]
        self.nodeDict = dict([(n,n) for n in self.nodeList])
        for i in range(grid.nHx):
            for j in range(grid.nHy):
                for k in range(grid.nHz):
                    #associate the element (i,j,k) with the
                    #left, front, bottom node
                    #get the left,front,bottom,node and its color
                    if (grid.getColor(i,j,k) == grid.black):
                        b0 = self.nodeList[grid.getNodeNumber(i,j,k)]
                        rx = self.nodeList[grid.getNodeNumber(i+1,j,k)]
                        ry = self.nodeList[grid.getNodeNumber(i,j+1,k)]
                        rz = self.nodeList[grid.getNodeNumber(i,j,k+1)]
                        r0 = self.nodeList[grid.getNodeNumber(i+1,j+1,k+1)]
                        bx = self.nodeList[grid.getNodeNumber(i,j+1,k+1)]
                        by = self.nodeList[grid.getNodeNumber(i+1,j,k+1)]
                        bz = self.nodeList[grid.getNodeNumber(i+1,j+1,k)]
                    else:
                        r0 = self.nodeList[grid.getNodeNumber(i,j,k)]
                        bx = self.nodeList[grid.getNodeNumber(i+1,j,k)]
                        by = self.nodeList[grid.getNodeNumber(i,j+1,k)]
                        bz = self.nodeList[grid.getNodeNumber(i,j,k+1)]
                        b0 = self.nodeList[grid.getNodeNumber(i+1,j+1,k+1)]
                        rx = self.nodeList[grid.getNodeNumber(i,j+1,k+1)]
                        ry = self.nodeList[grid.getNodeNumber(i+1,j,k+1)]
                        rz = self.nodeList[grid.getNodeNumber(i+1,j+1,k)]
                    self.newTetrahedron(nodes=[rx,by,bz,b0])
                    self.newTetrahedron(nodes=[ry,bz,bx,b0])
                    self.newTetrahedron(nodes=[rz,b0,bx,by])
                    self.newTetrahedron(nodes=[r0,bx,by,bz])
                    self.newTetrahedron(nodes=[b0,bx,by,bz])
        self.finalize()

    rectangularToTetrahedral = rectangularToTetrahedral6T

    def fixLocalNumbering(self):
        for TN in range(len(self.tetrahedronList)):
            self.tetrahedronList[TN].computeGeometricInfo()
            if edet(self.tetrahedronList[TN].linearMap) < 0:
                newNodes = list(self.tetrahedronList[TN].nodes)
                newNodes[2] = self.tetrahedronList[TN].nodes[1]
                newNodes[1] = self.tetrahedronList[TN].nodes[2]
                self.tetrahedronList[TN].nodes = newNodes
    def finalize(self):
        self.buildLists()
        #self.fixLocalNumbering()
        self.buildBoundaryMaps()
        self.buildArraysFromLists()
        self.hMax = 0.0
        self.hMin = 1.0e16
        self.sigmaMax = 0.0
        self.totalVolume = 0.0
        for T in self.tetrahedronList:
            T.computeGeometricInfo()
            self.hMax = max(T.diameter,self.hMax)
            self.hMin = min(T.diameter,self.hMin)
            self.sigmaMax = max(T.diameter/T.innerDiameter,self.sigmaMax)
            self.totalVolume += T.volume
    def buildLists(self):
        self.buildListsNodes()
        self.buildListsEdges()
        self.buildListsTriangles()
        self.buildListsTetrahedra()
        self.elementList = self.tetrahedronList
        self.elementBoundaryList = self.triangleList
    def buildListsNodes(self):
        keyList = self.nodeDict.keys()
        keyList.sort()
        self.nodeList=[]
        self.oldToNewNode=range(len(self.nodeDict))
        for nN,k in enumerate(keyList):
            self.oldToNewNode[self.nodeDict[k].N]=nN
            self.nodeDict[k].N = nN
            self.nodeList.append(self.nodeDict[k])

    def buildListsEdges(self):
        keyList = self.edgeDict.keys()
        keyList.sort()
        self.edgeList=[]
        for eN,k in enumerate(keyList):
            self.edgeDict[k].N = eN
            self.edgeList.append(self.edgeDict[k])

    def buildListsTriangles(self):
        keyList = self.triangleDict.keys()
        keyList.sort()
        self.triangleList=[]
        for tN,k in enumerate(keyList):
            self.triangleDict[k].N = tN
            self.triangleList.append(self.triangleDict[k])
        self.polygonList = self.triangleList

    def buildListsTetrahedra(self):
        keyList = self.tetrahedronDict.keys()
        keyList.sort()
        self.tetrahedronList=[]
        for TN,k in enumerate(keyList):
            self.tetrahedronDict[k].N = TN
            self.tetrahedronList.append(self.tetrahedronDict[k])
        self.polyhedronList = self.tetrahedronList

    def buildBoundaryMaps(self):
        """
        Extract a mapping tn -> list((TN,tnLocal)) that
        provides all elements with the boundary face (triangle) tn
        and the local triangle number for that triangle.
        Likewise build mappings for edges and nodes
        Also extract a list of the triangles with only one associate
        element; these are the external boundary triangles. Then extract
        the edges and nodes from the boundary triangles.
        """
        self.triangleMap=[[] for t in self.triangleList]
        self.edgeMap=[[] for e in self.edgeList]
        self.nodeMap=[[] for n in self.nodeList]
        self.boundaryTriangles=set()
        self.interiorTriangles=set()
        self.boundaryEdges=set()
        self.boundaryNodes=set()
        self.interiorEdges=set()
        self.interiorNodes=set()
        logEvent("Building triangle,edge, and node maps")
        for T in self.tetrahedronList:
            for localTriangleNumber,t in enumerate(T.triangles):
                self.triangleMap[t.N].append((T.N,localTriangleNumber))
            for localEdgeNumber,e in enumerate(T.edges):
                self.edgeMap[e.N].append((T.N,localEdgeNumber))
            for localNodeNumber,n in enumerate(T.nodes):
                self.nodeMap[n.N].append((T.N,localNodeNumber))
        logEvent("Extracting boundary and interior triangles")
        for tN,etList in enumerate(self.triangleMap):
            if len(etList) == 1:
                self.boundaryTriangles.add(self.triangleList[tN])
            else:
                self.interiorTriangles.add(self.triangleList[tN])
        logEvent("Extracting boundary edges and nodes")
        for t in self.boundaryTriangles:
            self.boundaryEdges.update(t.edges)
            self.boundaryNodes.update(t.nodes)
        logEvent("Extracting interior edges and nodes")
        for t in self.interiorTriangles:
            self.interiorEdges.update(t.edges)
            self.interiorNodes.update(t.nodes)
        self.boundaryMesh.buildFromSets(self.boundaryTriangles,
                                        self.boundaryEdges,self.boundaryNodes)

    def newTetrahedron(self,nodes):
        T = Tetrahedron(tetrahedronNumber=len(self.tetrahedronDict),
                        nodes=nodes)
        self.tetrahedronDict[T.nodes] = T
        self.registerTriangles(T)
        return T

    def registerEdges(self,t):
        for en,e in enumerate(t.edges):
            if self.edgeDict.has_key(e.nodes):
                t.edges[en]=self.edgeDict[e.nodes]
            else:
                eN=len(self.edgeDict)
                e.N=eN
                self.edgeDict[e.nodes]=e

    def registerTriangles(self,T):
        for tn,t in enumerate(T.triangles):
            if self.triangleDict.has_key(t.nodes):
                T.triangles[tn]=self.triangleDict[t.nodes]
            else:
                t.N=len(self.triangleDict)
                self.triangleDict[t.nodes]=t
                self.registerEdges(t)

    def registerNode(self,node):
        if self.nodeDict.has_key(node):
            node = self.nodeDict[node]
        else:
            node.N = len(self.nodeDict)
            self.nodeDict[node] = node
        return node

    def readMeshADH(self,filename,adhBase=1):
        meshIn = open(filename+'.3dm','r')
        firstLine = meshIn.readline()
        firstWords = firstLine.split()
        logEvent("Reading object=%s from file=%s" % (firstWords[0],filename))
        line = meshIn.readline()
        columns = line.split()
        tets = []
        tetEdges=set()
        tetTriangles=set()
        logEvent("Reading "+`filename`+" and building node lists for tetrahedra,triangles, and edges")
        #assume test are ordered by tet number
        while (columns[0] == 'E4T'):
            nodeNumbers = [int(c) - adhBase for c in columns[2:6]]
            nodeNumbers.sort()
            tets.append(array.array('i',nodeNumbers))
            tetTriangles.update([(nodeNumbers[1],nodeNumbers[2],nodeNumbers[3]),
                                 (nodeNumbers[0],nodeNumbers[2],nodeNumbers[3]),
                                 (nodeNumbers[0],nodeNumbers[1],nodeNumbers[3]),
                                 (nodeNumbers[0],nodeNumbers[1],nodeNumbers[2])])
            tetEdges.update([(nodeNumbers[0],nodeNumbers[1]),
                             (nodeNumbers[0],nodeNumbers[2]),
                             (nodeNumbers[0],nodeNumbers[3]),
                             (nodeNumbers[1],nodeNumbers[2]),
                             (nodeNumbers[1],nodeNumbers[3]),
                             (nodeNumbers[2],nodeNumbers[3])])
            line = meshIn.readline()
            columns = line.split()
        print "Building node list and dict"
        #assume nodes are ordered by node number
        while (len(columns) == 5):
            newNode = Node(int(columns[1]) - adhBase,
                           float(columns[2]),
                           float(columns[3]),
                           float(columns[4]))
            self.nodeList.append(newNode)
            self.nodeDict[newNode]=newNode
            line = meshIn.readline()
            columns = line.split()
        print "Number of tetrahedra:"+`len(tets)`
        print "Number of triangles :"+`len(tetTriangles)`
        print "Number of edges     :"+`len(tetEdges)`
        print "Number of nodes     :"+`len(self.nodeList)`
        print "Number of objects   :"+`len(tetEdges)+len(tetTriangles)+len(tets)+len(self.nodeList)`
        print "Building edge list"
        self.edgeList =[Edge(edgeNumber=eN,nodes=[self.nodeList[nN[0]],self.nodeList[nN[1]]]) \
                        for eN,nN in enumerate(tetEdges)]
        print "Building edge dict"
        self.edgeDict = dict([(e.nodes,e) for e in self.edgeList])
        print "Building triangle list"
        self.triangleList =[Triangle(triangleNumber=tN,nodes=[self.nodeList[nN[0]],self.nodeList[nN[1]],self.nodeList[nN[2]]],edgeDict=self.edgeDict) \
                            for tN,nN in enumerate(tetTriangles)]
        print "Building triangle dict"
        self.triangleDict = dict([(t.nodes,t) for t in self.triangleList])
        print "Building tetredron list"
        self.tetrahedronList = [Tetrahedron(tetrahedronNumber=TN,
                                            nodes=[self.nodeList[nN[0]],self.nodeList[nN[1]],self.nodeList[nN[2]],self.nodeList[nN[3]]],
                                            edgeDict=self.edgeDict,
                                            triangleDict=self.triangleDict) \
                                for TN,nN in enumerate(tets)]
        self.elementList = self.tetrahedronList
        self.elementBoundaryList = self.triangleList
        print "Building tetrahedron dict"
        self.tetrahedronDict = dict([(T.nodes,T) for T in self.tetrahedronList])
        print "Building boundary maps"
        self.buildBoundaryMaps()
    def writeMeshXdmf(self,ar,name='',t=0.0,init=False,meshChanged=False,tCount=0, EB=False):
        #print "Warning mwf hack for EB printing for tet writeMeshXdmf for now"
        #EB = True
        Mesh.writeMeshXdmf(self,ar,name,t,init,meshChanged,"Tetrahedron",tCount,EB=EB)
    def writeMeshEnsight(self,filename,description=None):
        base=1
        #write the casefile
        caseOut=open(filename+'.case','w')
        caseOut.write('FORMAT\n'+'type: ensight gold\n')
        caseOut.write('GEOMETRY\n'+'model: '+filename+'.geo\n')
        caseOut.close()
        meshOut=open(filename+'.geo','w')
        meshOut.write('Ensight Gold\n')
        meshOut.write('Unstructured Tetrahedral Mesh\n')
        meshOut.write('node id given\n')
        meshOut.write('element id given\n')
        #extents = 'extents\n %12.5E %12.5E\n %12.5E %12.5E\n %12.5E %12.5E\n' % (self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax)
        #meshOut.write('extents\n'+`self.xmin`+' '+`self.xmax`+'\n')
        meshOut.write('part \n'+'%10i\n' % 1)
        if description:
            meshOut.write(description+'\n')
        else:
            meshOut.write('A Mesh\n')
        meshOut.write('coordinates\n'+'%10i\n' % self.nNodes_global)
        for nN in range(self.nNodes_global):
            meshOut.write('%10i\n' % (nN+base))
        for nN in range(self.nNodes_global):
            meshOut.write('%12.5E\n' % self.nodeArray[nN,0])
        for nN in range(self.nNodes_global):
            meshOut.write('%12.5E\n' % self.nodeArray[nN,1])
        for nN in range(self.nNodes_global):
            meshOut.write('%12.5E\n' % self.nodeArray[nN,2])
        meshOut.write('tetra4\n'+'%10i\n' % self.nElements_global)
        for eN in range(self.nElements_global):
            meshOut.write('%10i\n' % (eN+base))
        for eN in range(self.nElements_global):
            meshOut.write('%10i%10i%10i%10i\n' % tuple((nN+base) for nN in self.elementNodesArray[eN,:]))
        meshOut.close()

    def appendMeshEnsight(self,filename,description=None):
        base=1
        #write the casefile
        caseOut=open(filename+'.case','a')
        caseOut.write('GEOMETRY\n'+'model: '+filename+'.geo\n')
        caseOut.close()
        meshOut=open(filename+'.geo','w')
        meshOut.write('Unstructured Tetrahedral Mesh\n\n')
        meshOut.write('node id given\n')
        meshOut.write('element id given\n')
        #extents = 'extents\n %12.5E %12.5E\n %12.5E %12.5E\n %12.5E %12.5E\n' % (self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax)
        #meshOut.write('extents\n'+`self.xmin`+' '+`self.xmax`+'\n')
        meshOut.write('part \n'+'%10i\n' % 1)
        if description:
            meshOut.write(description+'\n')
        else:
            meshOut.write("A Mesh")
        meshOut.write('coordinates\n'+'%10i\n' % len(self.nodeList))
        for n in self.nodeList:
            nN = n.N+base
            meshOut.write('%10i\n' % nN)
        for n in self.nodeList:
            meshOut.write('%12.5E\n' % n.p[X])
        for n in self.nodeList:
            meshOut.write('%12.5E\n' % n.p[Y])
        for n in self.nodeList:
            meshOut.write('%12.5E\n' % n.p[Z])
        meshOut.write('tetra4\n'+'%10i\n' % len(self.elementList))
        for e in self.elementList:
            eN = e.N + base
            meshOut.write('%10i\n' % eN)
        for e in self.elementList:
            meshOut.write('%10i%10i%10i%10i\n' % tuple(n.N+base for n in e.nodes))
        meshOut.close()

    def writeMeshADH(self,filename,adhBase=1):
        import cmeshTools
        cmeshTools.write3dmFiles(self.cmesh,filename,adhBase)

    def writeBoundaryFacesADH(self,filename,adhBase=1):
        boundaryFacesOut=open(filename,'w')
        for t in self.boundaryTriangles:
            TN = self.triangleMap[t.N][0][0]
            T = self.tetrahedronList[TN]
            localFaceNumber = self.triangleMap[t.N][0][1]
            T.computeGeometricInfo()
            DJ = edet(T.linearMap)
            if DJ < 0:
                #print "Negative determinant ="+`DJ`+" Swapping two nodes"
                newNodes = list(T.nodes)
                newNodes[3] = T.nodes[2]
                newNodes[2] = T.nodes[3]
                newBasis = [n - newNodes[0] for n in newNodes[1:]]
                newMap = ETen(newBasis[0],newBasis[1],newBasis[2])
                #print "New Determinant "+`edet(newMap)`
                if localFaceNumber == T.nodes[2]:
                    localFaceNumber = T.nodes[3]
                elif localFaceNumber == T.nodes[3]:
                    localFaceNumber = T.nodes[2]
            line = 'FCS %5i %5i %5i' % \
                          (T.N + adhBase,
                           localFaceNumber + adhBase,
                           1)
            #print line
            boundaryFacesOut.write(line+'\n')
        boundaryFacesOut.close()

    def writeBoundaryNodesADH(self,filename,adhBase=1):
        boundaryNodesOut=open(filename,'w')
        for n in self.boundaryNodes:
            line = 'NDS %5i %5i' % \
                          (n.N + adhBase,
                           1)
            #print line
            boundaryNodesOut.write(line+'\n')
        boundaryNodesOut.close()

    def refine4T(self,oldMesh):
        childrenDict={}
        for T in oldMesh.tetrahedronList:
            #deep copy old nodes because we'll renumber
            TNodes = [Node(eN,n.p[X],n.p[Y],n.p[Z]) for eN,n in enumerate(T.nodes)]
            for lnN,n in enumerate(TNodes): TNodes[lnN]=self.registerNode(n)
            #add new node
            T.computeGeometricInfo()
            newNode = Node(len(self.nodeDict),
                           T.barycenter[X],
                           T.barycenter[Y],
                           T.barycenter[Z])
            newNode = self.registerNode(newNode)
            T1=self.newTetrahedron([TNodes[0],TNodes[1],TNodes[2],newNode])
            T2=self.newTetrahedron([TNodes[1],TNodes[2],TNodes[3],newNode])
            T3=self.newTetrahedron([TNodes[2],TNodes[3],TNodes[0],newNode])
            T4=self.newTetrahedron([TNodes[3],TNodes[0],TNodes[1],newNode])
            childrenDict[T.N]=[T1,T2,T3,T4]
        self.finalize()
        return childrenDict

    def refineFreudenthalBey(self,oldMesh):
        logEvent("Refining the mesh using Freudenthal-Bey refinement")
        childrenDict={}
        for T in oldMesh.tetrahedronDict.values():
            #deep copy old nodes because we'll renumber
            TNodes = [Node(nN,n.p[X],n.p[Y],n.p[Z]) for nN,n in \
                      enumerate(T.nodes)]
            for lnN,n in enumerate(TNodes): TNodes[lnN]=self.registerNode(n)
            #add new nodes (midpoints of edges)
            #use local edge tuples as keys
            newNodes={}
            for et,en in T.edgeMap.iteritems():
                T.edges[en].computeGeometricInfo()
                p = T.edges[en].barycenter
                newNodes[et] = Node(en,p[X],p[Y],p[Z])

            #set the global node numbers
            for k,n in newNodes.iteritems(): newNodes[k]=self.registerNode(n)
            #add corner tets
            T1=self.newTetrahedron([TNodes[0],
                                 newNodes[(0,1)],
                                 newNodes[(0,2)],
                                 newNodes[(0,3)]])
            T2=self.newTetrahedron([TNodes[1],
                                 newNodes[(0,1)],
                                 newNodes[(1,2)],
                                 newNodes[(1,3)]])
            T3=self.newTetrahedron([TNodes[2],
                                 newNodes[(0,2)],
                                 newNodes[(1,2)],
                                 newNodes[(2,3)]])
            T4=self.newTetrahedron([TNodes[3],
                                 newNodes[(0,3)],
                                 newNodes[(1,3)],
                                 newNodes[(2,3)]])
            #add center tets
            #choose the shortest diagonal of the octahedron
            dLengths = [enorm(newNodes[(0,1)].p-newNodes[(2,3)].p),
                        enorm(newNodes[(0,2)].p-newNodes[(1,3)].p),
                        enorm(newNodes[(0,3)].p-newNodes[(1,2)].p)]
            shortestEdgeLength = min(dLengths)
            if shortestEdgeLength == dLengths[0]:
                #diagonal (0,1)(2,3)
                T5=self.newTetrahedron([newNodes[(0,1)],
                                        newNodes[(2,3)],
                                        newNodes[(0,3)],
                                        newNodes[(1,3)]])
                T6=self.newTetrahedron([newNodes[(0,1)],
                                        newNodes[(2,3)],
                                        newNodes[(0,3)],
                                        newNodes[(0,2)]])
                T7=self.newTetrahedron([newNodes[(0,1)],
                                        newNodes[(2,3)],
                                        newNodes[(0,2)],
                                        newNodes[(1,2)]])
                T8=self.newTetrahedron([newNodes[(0,1)],
                                        newNodes[(2,3)],
                                        newNodes[(1,2)],
                                        newNodes[(1,3)]])
            elif shortestEdgeLength == dLengths[1]:
                #diagonal (0,2)(1,3)
                T5=self.newTetrahedron([newNodes[(0,2)],
                                        newNodes[(1,3)],
                                        newNodes[(0,3)],
                                        newNodes[(2,3)]])
                T6=self.newTetrahedron([newNodes[(0,2)],
                                        newNodes[(1,3)],
                                        newNodes[(2,3)],
                                        newNodes[(1,2)]])
                T7=self.newTetrahedron([newNodes[(0,2)],
                                        newNodes[(1,3)],
                                        newNodes[(1,2)],
                                        newNodes[(0,1)]])
                T8=self.newTetrahedron([newNodes[(0,2)],
                                        newNodes[(1,3)],
                                        newNodes[(0,1)],
                                        newNodes[(0,3)]])
            else:
                #diagonal (0,3)(1,2)
                T5=self.newTetrahedron([newNodes[(0,3)],
                                        newNodes[(1,2)],
                                        newNodes[(0,1)],
                                        newNodes[(1,3)]])
                T6=self.newTetrahedron([newNodes[(0,3)],
                                        newNodes[(1,2)],
                                        newNodes[(1,3)],
                                        newNodes[(2,3)]])
                T7=self.newTetrahedron([newNodes[(0,3)],
                                        newNodes[(1,2)],
                                        newNodes[(2,3)],
                                        newNodes[(0,2)]])
                T8=self.newTetrahedron([newNodes[(0,3)],
                                        newNodes[(1,2)],
                                        newNodes[(0,2)],
                                        newNodes[(0,1)]])
            childrenDict[T.N]=[T1,T2,T3,T4,T5,T6,T7,T8]
        self.finalize()
        return childrenDict
        #for debugging: print each tet
        #self.edgeList=[]
        #Tlist = self.tetrahedronDict.values()
        #for T in Tlist:
        #    self.edgeList = self.edgeList + T.edges

    def refine(self,oldMesh):
        return self.refineFreudenthalBey(oldMesh)

    def generateFromTetgenFiles(self,filebase,base,skipGeometricInit=True,parallel=False):
        import cmeshTools
        logEvent(memory("declaring CMesh"),level=4)
        self.cmesh = cmeshTools.CMesh()
        logEvent(memory("Initializing CMesh"),level=4)
        if parallel:
            cmeshTools.generateFromTetgenFilesParallel(self.cmesh,filebase,base)
        else:
            cmeshTools.generateFromTetgenFiles(self.cmesh,filebase,base)
        logEvent(memory("calling cmeshTools.generateFromTetgenFiles","cmeshTools"),level=4)
        if skipGeometricInit == False:
            cmeshTools.allocateGeometricInfo_tetrahedron(self.cmesh)
            cmeshTools.computeGeometricInfo_tetrahedron(self.cmesh)
        self.buildFromC(self.cmesh)
        logEvent(memory("calling buildFromC"),level=4)
    def generateFrom3DMFile(self,filebase,base=1):
        import cmeshTools
        self.cmesh = cmeshTools.CMesh()
        cmeshTools.generateFrom3DMFile(self.cmesh,filebase,base)
        cmeshTools.allocateGeometricInfo_tetrahedron(self.cmesh)
        cmeshTools.computeGeometricInfo_tetrahedron(self.cmesh)
        self.buildFromC(self.cmesh)
    def writeTetgenFiles(self,filebase,base):
        import cmeshTools
        cmeshTools.writeTetgenFiles(self.cmesh,filebase,base)
    def meshInfo(self):
        minfo = """Number of tetrahedra : %d
Number of triangles  : %d
Number of edges      : %d
Number of nodes      : %d
max(sigma_k)         : %f
min(h_k)             : %f\n""" % (self.nElements_global,
                                  self.nElementBoundaries_global,
                                  self.nEdges_global,
                                  self.nNodes_global,
                                  self.sigmaMax,
                                  self.hMin)
        if self.subdomainMesh != self:
            sinfo = self.subdomainMesh.meshInfo()
            info = "*** Global ***\n" + minfo + "\n*** Local ***\n" + sinfo
            return info
        return minfo

class HexahedralMesh(Mesh):
    """A mesh of hexahedra.

    """

    def __init__(self):
        Mesh.__init__(self)
        self.nodeDict={}
        self.edgeDict={}
        self.faceDict={}
        self.faceList=[]
        self.elemDict={}
        self.elemList=[]
        self.oldToNewNode=[]
        self.boundaryMesh=QuadrilateralMesh()
    def meshType(self):
        return 'cuboid'
    def computeGeometricInfo(self):
        import cmeshTools
        print "no info yet for hexahedral mesh"
        #cmeshTools.computeGeometricInfo_tetrahedron(self.cmesh)
    def generateHexahedralMeshFromRectangularGrid(self,nx,ny,nz,Lx,Ly,Lz):
        import cmeshTools
        self.cmesh = cmeshTools.CMesh()
        cmeshTools.generateHexahedralMeshFromRectangularGrid(nx,ny,nz,0,0,0,Lx,Ly,Lz,self.cmesh)
        cmeshTools.allocateGeometricInfo_hexahedron(self.cmesh)
        cmeshTools.computeGeometricInfo_hexahedron(self.cmesh)
        self.buildFromC(self.cmesh)

    def finalize(self):
        self.buildLists()
        #self.fixLocalNumbering()
        self.buildBoundaryMaps()
        self.buildArraysFromLists()
        self.hMax = 0.0
        self.hMin = 1.0e16
        self.sigmaMax = 0.0
        self.totalVolume = 0.0
        for T in self.tetrahedronList:
            T.computeGeometricInfo()
            self.hMax = max(T.diameter,self.hMax)
            self.hMin = min(T.diameter,self.hMin)
            self.sigmaMax = max(T.diameter/T.innerDiameter,self.sigmaMax)
            self.totalVolume += T.volume

    def buildLists(self):
        self.buildListsNodes()
        self.buildListsEdges()
        self.buildListsFaces()
        self.buildListsElems()
        self.elementList = self.elemList
        self.elementBoundaryList = self.faceList
    def buildListsNodes(self):
        keyList = self.nodeDict.keys()
        keyList.sort()
        self.nodeList=[]
        self.oldToNewNode=range(len(self.nodeDict))
        for nN,k in enumerate(keyList):
            self.oldToNewNode[self.nodeDict[k].N]=nN
            self.nodeDict[k].N = nN
            self.nodeList.append(self.nodeDict[k])

    def buildListsEdges(self):
        keyList = self.edgeDict.keys()
        keyList.sort()
        self.edgeList=[]
        for eN,k in enumerate(keyList):
            self.edgeDict[k].N = eN
            self.edgeList.append(self.edgeDict[k])

    def buildListsFaces(self):
        keyList = self.faceDict.keys()
        keyList.sort()
        self.triangleList=[]
        for tN,k in enumerate(keyList):
            self.faceDict[k].N = tN
            self.faceList.append(self.faceDict[k])
        self.polygonList = self.faceList

    def buildListsElems(self):
        keyList = self.elemDict.keys()
        keyList.sort()
        self.elemList=[]
        for TN,k in enumerate(keyList):
            self.elemDict[k].N = TN
            self.elemList.append(self.elemDict[k])
        self.polyhedronList = self.elemList

    def buildBoundaryMaps(self):
        """
        Extract a mapping tn -> list((TN,tnLocal)) that
        provides all elements with the boundary face  tn
        and the local triangle number for that face
        Likewise build mappings for edges and nodes
        Also extract a list of the triangles with only one associate
        element; these are the external boundary triangles. Then extract
        the edges and nodes from the boundary triangles.
        """
        self.faceMap=[[] for t in self.faceList]
        self.edgeMap=[[] for e in self.edgeList]
        self.nodeMap=[[] for n in self.nodeList]
        self.boundaryTriangles=set()
        self.interiorTriangles=set()
        self.boundaryEdges=set()
        self.boundaryNodes=set()
        self.interiorEdges=set()
        self.interiorNodes=set()
        logEvent("Building triangle,edge, and node maps")
        for T in self.elemList:
            for localFaceNumber,t in enumerate(T.faces):
                self.faceMap[t.N].append((T.N,localFaceNumber))
            for localEdgeNumber,e in enumerate(T.edges):
                self.edgeMap[e.N].append((T.N,localEdgeNumber))
            for localNodeNumber,n in enumerate(T.nodes):
                self.nodeMap[n.N].append((T.N,localNodeNumber))
        logEvent("Extracting boundary and interior triangles")
        for tN,etList in enumerate(self.faceMap):
            if len(etList) == 1:
                self.boundaryFaces.add(self.faceList[tN])
            else:
                self.interiorFaces.add(self.faceList[tN])
        logEvent("Extracting boundary edges and nodes")
        for t in self.boundaryTriangles:
            self.boundaryEdges.update(t.edges)
            self.boundaryNodes.update(t.nodes)
        logEvent("Extracting interior edges and nodes")
        for t in self.interiorTriangles:
            self.interiorEdges.update(t.edges)
            self.interiorNodes.update(t.nodes)
        self.boundaryMesh.buildFromSets(self.boundaryFaces,
                                        self.boundaryEdges,self.boundaryNodes)

    def registerEdges(self,t):
        for en,e in enumerate(t.edges):
            if self.edgeDict.has_key(e.nodes):
                t.edges[en]=self.edgeDict[e.nodes]
            else:
                eN=len(self.edgeDict)
                e.N=eN
                self.edgeDict[e.nodes]=e

    def registerFaces(self,T):
        for tn,t in enumerate(T.faces):
            if self.faceDict.has_key(t.nodes):
                T.faces[tn]=self.faceDict[t.nodes]
            else:
                t.N=len(self.faceDict)
                self.faceDict[t.nodes]=t
                self.registerEdges(t)

    def registerNode(self,node):
        if self.nodeDict.has_key(node):
            node = self.nodeDict[node]
        else:
            node.N = len(self.nodeDict)
            self.nodeDict[node] = node
        return node

#    def refine(self,oldMesh):
#        return self.refineFreudenthalBey(oldMesh)

    def meshInfo(self):
        minfo = """Number of hexahedra  : %d
Number of faces      : %d
Number of edges      : %d
Number of nodes      : %d
max(sigma_k)         : %d
min(h_k)             : %d\n""" % (self.nElements_global,
                                  self.nElementBoundaries_global,
                                  self.nEdges_global,
                                  self.nNodes_global,
                                  self.sigmaMax,
                                  self.hMin)
        if self.subdomainMesh != self:
            sinfo = self.subdomainMesh.meshInfo()
            info = "*** Global ***\n" + minfo + "\n*** Local ***\n" + sinfo
            return info
        return minfo

    def writeMeshXdmf(self,ar,name='',t=0.0,init=False,meshChanged=False,tCount=0,EB=False):
        Mesh.writeMeshXdmf(self,ar,name,t,init,meshChanged,"Hexahedron",tCount,EB=EB)

    def generateFromHexFile(self,filebase,base=0):
        import cmeshTools
        self.cmesh = cmeshTools.CMesh()
        cmeshTools.generateFromHexFile(self.cmesh,filebase,base)
        cmeshTools.allocateGeometricInfo_hexahedron(self.cmesh)
        cmeshTools.computeGeometricInfo_hexahedron(self.cmesh)
        self.buildFromC(self.cmesh)

class Mesh2DM(Mesh):
    """A triangular mesh based on an ADH 3dm file"""
    def __init__(self,filename,adhBase=1):
        meshIn = open(filename+'.3dm','r')
        firstLine = meshIn.readline()
        firstWords = firstLine.split()
        logEvent("Reading object=%s from file=%s" % (firstWords[0],filename))
        line = meshIn.readline()
        columns = line.split()
        #read in the tetrahedra and nodes as memory-efficiently as possible
        tn0 = array.array('i')
        tn1 = array.array('i')
        tn2 = array.array('i')
        material = array.array('i')
        nx  = array.array('d')
        ny  = array.array('d')
        nz  = array.array('d')
        print "Reading "+`filename`
        #assume tets are ordered by tet number
        while (len(columns) > 0 and (columns[0] == 'E3T' or columns[0] == 'GE3')):
            tn0.append(int(columns[2]))
            tn1.append(int(columns[3]))
            tn2.append(int(columns[4]))
            material.append(int(columns[5]))
            line = meshIn.readline()
            columns = line.split()
        #allow for missing lines
        while (len(columns) == 0):
            line = meshIn.readline()
            columns = line.split()
        #assume nodes are ordered by node number
        while (len(columns) == 5):
            nx.append(float(columns[2]))
            ny.append(float(columns[3]))
            nz.append(float(columns[4]))
            line = meshIn.readline()
            columns = line.split()
        meshIn.close()
        print "Allocating node and element arrays"
        self.nTriangles_global = len(tn0)
        self.triangleArray = np.zeros(
            (self.nTriangles_global,3),'i')
        tA = self.triangleArray
        self.triangleMaterialArray = np.zeros(
            (self.nTriangles_global,),'i')
        tMA = self.triangleMaterialArray
        self.nNodes_global = len(nx)
        self.nodeArray = np.zeros((self.nNodes_global,3),'d')
        for tN in range(self.nTriangles_global):
            tA[tN,0] = tn0[tN] - adhBase
            tA[tN,1] = tn1[tN] - adhBase
            tA[tN,2] = tn2[tN] - adhBase
            tMA[tN]  = material[tN] - adhBase
        for nN in range(self.nNodes_global):
            self.nodeArray[nN,0]= nx[nN]
            self.nodeArray[nN,1]= ny[nN]
            self.nodeArray[nN,2]= nz[nN]
        print "Deleting temporary storage"
        del tn0,tn1,tn2,nx,ny,nz
        self.nElements_global = self.nTriangles_global
        self.elementNodesArray = self.triangleArray
        self.elementMaterialTypes = self.triangleMaterialArray
        print "Number of triangles:"+`self.nElements_global`
        print "Number of nodes     :"+`self.nNodes_global`
        #archive with Xdmf
        self.nNodes_element = 3
        self.arGridCollection = None
        self.arGrid = None; self.arTime = None

    def buildEdgeArrays(self):
        print "Extracting edges triangles dictionary"
        edges_triangles={}
        t=self.triangleArray
        self.nInteriorEdges_global=0
        for N in range(self.nTriangles_global):
            #sort node numbers so the nodes can
            #uniquely identify the triangles/edges
            n = list(t[N,:])
            n.sort()
            edges = [(n[0],n[1]),
                     (n[0],n[2]),
                     (n[1],n[2])]
            for t in triangles:
                if edges_triangles.has_key(t):
                    edges_triangles[t].append(N)
                    self.nInteriorTriangles_global+=1
                else:
                    edges_triangles[t]=[N]
        print "Building edge and exterior arrays"
        self.nEdges_global = len(edges_triangles)
        self.edgeArray = np.zeros(
            (self.nEdges_global,2),'i')
        self.edgeMaterialArray = np.zeros(
            (self.nEdges_global,2),'i')
        self.interiorEdgeArray = np.zeros(
            (self.nInteriorEdges_global,),'i')
        self.nExteriorEdges_global = self.nEdges_global - \
                                     self.nInteriorEdges_global
        self.exteriorEdgeArray = np.zeros(
            (self.nExteriorEdges_global,),'i')
        eN=0
        ieN=0
        eeN=0
        exteriorNodes=set()
        eA = self.edgeArray
        eMA = self.edgeMaterialArray
        tMA = self.triangleMaterialArray
        for eNodes,tlist in edges_triangles.iteritems():
            eA[eN,0]=eNodes[0]
            eA[eN,1]=eNodes[1]
            if len(tlist)==2:
                self.interiorEdgeArray[ieN]=eN
                eMA[eN][0]= tMA[tlist[0]]
                eMA[eN][1]= tMA[Tlist[1]]
                ieN+=1
            else:
                exteriorNodes.update(tNodes)
                self.exteriorEdgeArray[eeN]=eN
                eMA[eN][0]=tMA[tlist[0]]
                eeN+=1
            eN+=1
        self.nExteriorNodes_global = len(exteriorNodes)
        self.exteriorNodeArray = np.zeros(
            (self.nExteriorNodes_global,),'i')
        self.globalToExteriorNodeArray = np.zeros(
            (self.nNodes_global,),'i')
        for nExtN,nN in enumerate(exteriorNodes):
            self.exteriorNodeArray[nExtN]=nN
            self.globalToExteriorNodeArray[nN]=nExtN
        print "Number of edges         :"+`self.nEdges_global`
        print "Number on interior      :"+`self.nInteriorEdges_global`
        print "Number on exterior      :"+`self.nExteriorEdges_global`
        print "Number of exterior nodes:"+`self.nExteriorNodes_global`
        #at this point we can easily build a boundary mesh by renumbering using
        #exteriorNodeArray and exteriorEdgeArray to renumber
        #and the info in nodeArray and edgeArray

    def writeBoundaryMeshADH(self,filename,adhBase=1):
        #I'll print it using node numbers from the 3D mesh
        meshOut = open(filename+'Boundary.3dm','w')
        meshOut.write('MESH1D\n')
        for eeN in range(self.nExteriorEdges_global):
            eN = self.exteriorEdgeArray[eeN]
            n0 = self.edgeArray[eN][0] + adhBase
            n1 = self.edgeArray[eN][1] + adhBase
            m  = self.edgeMaterialArray[eN][0] + adhBase
            line = 'E3T %5i %5i %5i %5i' % \
                          (tN+adhBase,n0,n1,m)
            meshOut.write(line+'\n')
        meshOut.close()

    def writeMeshEnsight(self,filename,description=None):
        base=1
        #write the casefile
        caseOut=open(filename+'.case','w')
        caseOut.write('FORMAT\n'+'type: ensight gold\n')
        caseOut.write('GEOMETRY\n'+'model: '+filename+'.geo\n')
        caseOut.close()
        meshOut=open(filename+'.geo','w')
        meshOut.write('Ensight Gold\n')
        meshOut.write('Unstructured Triangular Mesh\n')
        meshOut.write('node id given\n')
        meshOut.write('element id given\n')
        #extents = 'extents\n %12.5E %12.5E\n %12.5E %12.5E\n %12.5E %12.5E\n' % (self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax)
        #meshOut.write('extents\n'+`self.xmin`+' '+`self.xmax`+'\n')
        meshOut.write('part \n'+'%10i\n' % 1)
        if description:
            meshOut.write(description+'\n')
        else:
            meshOut.write('A Mesh\n')
        meshOut.write('coordinates\n'+'%10i\n' % self.nNodes_global)
        for nN in range(self.nNodes_global):
            ensightNodeNumber = (nN+base)
            meshOut.write('%10i\n' % ensightNodeNumber)
        for nN in range(self.nNodes_global):
            meshOut.write('%12.5E\n' % self.nodeArray[nN,0])
        for nN in range(self.nNodes_global):
            meshOut.write('%12.5E\n' % self.nodeArray[nN,1])
        for nN in range(self.nNodes_global):
            meshOut.write('%12.5E\n' % self.nodeArray[nN,2])
        meshOut.write('tria3\n'+'%10i\n' % self.nTriangles_global)
        for tN in range(self.nTriangles_global):
            ensightElementNumber = tN + base
            meshOut.write('%10i\n' % ensightElementNumber)
        tA = self.triangleArray
        for tN in range(self.nTriangles_global):
            meshOut.write('%10i%10i%10i\n' % (tA[tN,0]+base,tA[tN,1]+base,tA[tN,2]+base))
        meshOut.close()

    def writeMeshXdmf(self,ar,name='',t=0.0,init=False,meshChanged=False,Xdmf_ElementTopology="Triangle",tCount=0):
        if self.arGridCollection is not None:
            init = False
        elif not init:
            self.arGridCollection = ar.domain.find("Grid")
        if init:
            self.arGridCollection = SubElement(ar.domain,"Grid",{"Name":"Mesh "+name,
                                                                "GridType":"Collection",
                                                                "CollectionType":"Temporal"})
        if self.arGrid is None or self.arTime.get('Value') != str(t):
            #
            #topology and geometry
            #
            if ar.global_sync:
                self.arGrid = SubElement(self.arGridCollection,"Grid",{"GridType":"Uniform"})
                self.arTime = SubElement(self.arGrid,"Time",{"Value":str(t),"Name":str(tCount)})
                topology = SubElement(self.arGrid,"Topology",
                                      {"Type":Xdmf_ElementTopology,
                                       "NumberOfElements":str(self.globalMesh.nElements_global)})
                elements = SubElement(topology,"DataItem",
                                      {"Format":ar.dataItemFormat,
                                       "DataType":"Int",
                                       "Dimensions":"%i %i" % (self.globalMesh.nElements_global,
                                                               self.nNodes_element)})
                geometry = SubElement(self.arGrid,"Geometry",{"Type":"XYZ"})
                nodes    = SubElement(geometry,"DataItem",
                                      {"Format":ar.dataItemFormat,
                                       "DataType":"Float",
                                       "Precision":"8",
                                       "Dimensions":"%i %i" % (self.globalMesh.nNodes_global,3)})
                #material types
                elementMaterialTypes = SubElement(self.arGrid,"Attribute",{"Name":"elementMaterialTypes",
                                                                           "AttributeType":"Scalar",
                                                                           "Center":"Cell"})
                elementMaterialTypesValues = SubElement(elementMaterialTypes,"DataItem",
                                                        {"Format":ar.dataItemFormat,
                                                         "DataType":"Int",
                                                         "Dimensions":"%i" % (self.globalMesh.nElements_global,)})
                if ar.hdfFile is not None:
                    if ar.has_h5py:
                        elements.text = ar.hdfFilename+":/elements"+name+`tCount`
                        nodes.text = ar.hdfFilename+":/nodes"+name+`tCount`
                        elementMaterialTypesValues.text = ar.hdfFilename+":/"+"elementMaterialTypes"+"_t"+str(tCount)
                        if init or meshChanged:
                            ar.create_dataset_sync('elements'+name+`tCount`,
                                                   offsets = self.globalMesh.elementOffsets_subdomain_owned,
                                                   data = self.globalMesh.nodeNumbering_subdomain2global[self.elementNodesArray[:self.nElements_owned]])
                            ar.create_dataset_sync('nodes'+name+`tCount`,
                                                   offsets = self.globalMesh.nodeOffsets_subdomain_owned,
                                                   data = self.nodeArray[:self.nNodes_owned])
                            ar.create_dataset_sync("elementMaterialTypes"+"_t"+str(tCount),
                                                   offsets = self.globalMesh.elementOffsets_subdomain_owned,
                                                   data = self.elementMaterialTypes[:self.nElements_owned])
                    else:
                        assert False, "global_sync with pytables not supported"
                else:
                    assert False, "global_sync with text heavy data not supported"
            else:
                self.arGrid = SubElement(self.arGridCollection,"Grid",{"GridType":"Uniform"})
                self.arTime = SubElement(self.arGrid,"Time",{"Value":str(t),"Name":str(tCount)})
                topology = SubElement(self.arGrid,"Topology",
                                      {"Type":Xdmf_ElementTopology,
                                       "NumberOfElements":str(self.nElements_owned)})
                elements = SubElement(topology,"DataItem",
                                      {"Format":ar.dataItemFormat,
                                       "DataType":"Int",
                                       "Dimensions":"%i %i" % (self.nElements_owned,self.nNodes_element)})
                geometry = SubElement(self.arGrid,"Geometry",{"Type":"XYZ"})
                nodes    = SubElement(geometry,"DataItem",
                                      {"Format":ar.dataItemFormat,
                                       "DataType":"Float",
                                       "Precision":"8",
                                       "Dimensions":"%i %i" % (self.nNodes_global,3)})
                #material types
                elementMaterialTypes = SubElement(self.arGrid,"Attribute",{"Name":"elementMaterialTypes",
                                                                           "AttributeType":"Scalar",
                                                                           "Center":"Cell"})
                elementMaterialTypesValues = SubElement(elementMaterialTypes,"DataItem",
                                                        {"Format":ar.dataItemFormat,
                                                         "DataType":"Int",
                                                         "Dimensions":"%i" % (self.nElements_owned,)})
                if ar.hdfFile is not None:
                    if ar.has_h5py:
                        elements.text = ar.hdfFilename+":/elements"+`ar.comm.rank()`+name+`tCount`
                        nodes.text = ar.hdfFilename+":/nodes"+`ar.comm.rank()`+name+`tCount`
                        elementMaterialTypesValues.text = ar.hdfFilename+":/"+"elementMaterialTypes"+"_p"+`ar.comm.rank()`+"_t"+str(tCount)
                        if init or meshChanged:
                            ar.create_dataset_async('elements'+`ar.comm.rank()`+name+`tCount`, data = self.elementNodesArray[:self.nElements_owned])
                            ar.create_dataset_async('nodes'+`ar.comm.rank()`+name+`tCount`, data = self.nodeArray)
                            ar.create_dataset_async("elementMaterialTypes"+"_p"+`ar.comm.rank()`+"_t"+str(tCount), data = self.elementMaterialTypes[:self.nElements_owned])
                    else:
                        elements.text = ar.hdfFilename+":/elements"+name+`tCount`
                        nodes.text = ar.hdfFilename+":/nodes"+name+`tCount`
                        elementMaterialTypesValues.text = ar.hdfFilename+":/"+"elementMaterialTypes"+str(tCount)
                        if init or meshChanged:
                            ar.hdfFile.createArray("/",'elements'+name+`tCount`,self.elementNodesArray[:self.nElements_owned])
                            ar.hdfFile.createArray("/",'nodes'+name+`tCount`,self.nodeArray)
                            ar.hdfFile.createArray("/","elementMaterialTypes"+str(tCount),self.elementMaterialTypes[:self.nElements_owned])
                else:
                    SubElement(elements,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/elements"+name+".txt"})
                    SubElement(nodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodes"+name+".txt"})
                    SubElement(elementMaterialTypesValues,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+"elementMaterialTypes"+str(tCount)+".txt"})
                    if init or meshChanged:
                        np.savetxt(ar.textDataDir+"/elements"+name+".txt",self.elementNodesArray[:self.nElements_owned],fmt='%d')
                        np.savetxt(ar.textDataDir+"/nodes"+name+".txt",self.nodeArray)
                        np.savetxt(ar.textDataDir+"/"+"elementMaterialTypes"+str(tCount)+".txt",self.elementMaterialTypes[:self.nElements_owned])


class Mesh3DM(Mesh):
    """
    A Mesh for reading in tetrahedral meshes in the .3dm format
    """
    def __init__(self,filename,adhBase=1):
        meshIn = open(filename+'.3dm','r')
        firstLine = meshIn.readline()
        firstWords = firstLine.split()
        print "Reading object=%s from file=%s" % (firstWords[0],filename)
        line = meshIn.readline()
        columns = line.split()
        #read in the tetrahedra and nodes as memory-efficiently as possible
        Tn0 = array.array('i')
        Tn1 = array.array('i')
        Tn2 = array.array('i')
        Tn3 = array.array('i')
        material = array.array('i')
        nx  = array.array('d')
        ny  = array.array('d')
        nz  = array.array('d')
        print "Reading "+`filename`
        #assume tets are ordered by tet number
        while (len(columns) > 0 and (columns[0] == 'E4T' or columns[0] == 'GE4')):
            Tn0.append(int(columns[2]))
            Tn1.append(int(columns[3]))
            Tn2.append(int(columns[4]))
            Tn3.append(int(columns[5]))
            material.append(int(columns[6]))
            line = meshIn.readline()
            columns = line.split()
        #assume nodes are ordered by node number
        while (len(columns) == 5):
            nx.append(float(columns[2]))
            ny.append(float(columns[3]))
            nz.append(float(columns[4]))
            line = meshIn.readline()
            columns = line.split()
        meshIn.close()
        print "Allocating node and element arrays"
        self.nTetrahedra_global = len(Tn0)
        self.tetrahedronArray = np.zeros(
            (self.nTetrahedra_global,4),'i')
        TA = self.tetrahedronArray
        self.tetrahedronMaterialArray = np.zeros(
            (self.nTetrahedra_global,),'i')
        TMA = self.tetrahedronMaterialArray
        self.nNodes_global = len(nx)
        self.nodeArray = np.zeros((self.nNodes_global,3),'d')
        for TN in range(self.nTetrahedra_global):
            TA[TN,0] = Tn0[TN] - adhBase
            TA[TN,1] = Tn1[TN] - adhBase
            TA[TN,2] = Tn2[TN] - adhBase
            TA[TN,3] = Tn3[TN] - adhBase
            TMA[TN]  = material[TN] - adhBase
        for nN in range(self.nNodes_global):
            self.nodeArray[nN,0]= nx[nN]
            self.nodeArray[nN,1]= ny[nN]
            self.nodeArray[nN,2]= nz[nN]
        print "Deleting temporary storage"
        del Tn0,Tn1,Tn2,Tn3,nx,ny,nz
        self.nElements_global = self.nTetrahedra_global
        self.elementNodesArray = self.tetrahedronArray
        self.elementMaterialTypes = self.tetrahedronMaterialArray
        self.arGridCollection=None
        print "Number of tetrahedra:"+`self.nElements_global`
        print "Number of nodes     :"+`self.nNodes_global`

    def buildTriangleArrays(self):
        print "Extracting triangles tetrahedra dictionary"
        triangles_tetrahedra={}
        T=self.tetrahedronArray
        self.nInteriorTriangles_global=0
        for N in range(self.nTetrahedra_global):
            #sort node numbers so the nodes can
            #uniquely identify the triangles/edges
            n = list(T[N,:])
            n.sort()
            triangles = [(n[0],n[1],n[2]),
                         (n[0],n[1],n[3]),
                         (n[0],n[2],n[3]),
                         (n[1],n[2],n[3])]
            for t in triangles:
                if triangles_tetrahedra.has_key(t):
                    triangles_tetrahedra[t].append(N)
                    self.nInteriorTriangles_global+=1
                else:
                    triangles_tetrahedra[t]=[N]
        print "Building triangle and exterior arrays"
        self.nTriangles_global = len(triangles_tetrahedra)
        self.triangleArray = np.zeros(
            (self.nTriangles_global,3),'i')
        self.triangleMaterialArray = np.zeros(
            (self.nTriangles_global,2),'i')
        self.interiorTriangleArray = np.zeros(
            (self.nInteriorTriangles_global,),'i')
        self.nExteriorTriangles_global = self.nTriangles_global - \
                                         self.nInteriorTriangles_global
        self.exteriorTriangleArray = np.zeros(
            (self.nExteriorTriangles_global,),'i')
        tN=0
        itN=0
        etN=0
        exteriorNodes=set()
        tA = self.triangleArray
        tMA = self.triangleMaterialArray
        TMA = self.tetrahedronMaterialArray
        for tNodes,Tlist in triangles_tetrahedra.iteritems():
            tA[tN,0]=tNodes[0]
            tA[tN,1]=tNodes[1]
            tA[tN,2]=tNodes[2]
            if len(Tlist)==2:
                self.interiorTriangleArray[itN]=tN
                tMA[tN][0]= TMA[Tlist[0]]
                tMA[tN][1]= TMA[Tlist[1]]
                itN+=1
            else:
                exteriorNodes.update(tNodes)
                self.exteriorTriangleArray[etN]=tN
                tMA[tN][0]=TMA[Tlist[0]]
                etN+=1
            tN+=1
        self.nExteriorNodes_global = len(exteriorNodes)
        self.exteriorNodeArray = np.zeros(
            (self.nExteriorNodes_global,),'i')
        self.globalToExteriorNodeArray = np.zeros(
            (self.nNodes_global,),'i')
        for nExtN,nN in enumerate(exteriorNodes):
            self.exteriorNodeArray[nExtN]=nN
            self.globalToExteriorNodeArray[nN]=nExtN
        print "Number of triangles     :"+`self.nTriangles_global`
        print "Number on interior      :"+`self.nInteriorTriangles_global`
        print "Number on exterior      :"+`self.nExteriorTriangles_global`
        print "Number of exterior nodes:"+`self.nExteriorNodes_global`
        #at this point we can easily build a boundary mesh by renumbering using
        #exteriorNodeArray and exteriorTriangleArray to renumber
        #and the info in nodeArray and triangleArray

    def buildEdgeArray(self):
        print "Extracting set of edges"
        edges = set()
        t=self.triangleArray
        for N in range(self.nTriangles_global):
            #triangle nodes are assumed sorted
            edges.update([(t[N,0],t[N,1]),
                          (t[N,0],t[N,2]),
                          (t[N,1],t[N,2])])
        print "Building edgeArray"
        self.nEdges_global = len(edges)
        self.edgeArray = np.zeros(
            (self.nEdges_global,2),'i')
        eN=0
        for e in edges:
            self.edgeArray[eN][0] = e[0]
            self.edgeArray[eN][1] = e[1]
        del edges
        print "Number of edges     :"+`self.nEdges_global`

    def writeBoundaryMeshADH(self,filename,adhBase=1):
        #I'll print it using node numbers from the 3D mesh
        meshOut = open(filename+'Boundary.3dm','w')
        meshOut.write('MESH2D\n')
        for tN in self.exteriorTriangleArray:
            n0 = self.triangleArray[tN][0] + adhBase
            n1 = self.triangleArray[tN][1] + adhBase
            n2 = self.triangleArray[tN][2] + adhBase
            m  = self.triangleMaterialArray[tN][0] + adhBase
            line = 'E3T %5i %5i %5i %5i %5i' % \
                          (tN+adhBase,n0,n1,n2,m)
            meshOut.write(line+'\n')
        for nN in self.exteriorNodeArray:
            n  = self.nodeArray[nN]
            line = 'ND %5i %14.8e %14.8e %14.8e' % \
                          (nN + adhBase,n[0],n[1],n[2])
            #print line
            meshOut.write(line+'\n')
        meshOut.close()

    def writeMeshEnsight(self,filename,description=None):
        base=1
        #write the casefile
        caseOut=open(filename+'.case','w')
        caseOut.write('FORMAT\n'+'type: ensight gold\n')
        caseOut.write('GEOMETRY\n'+'model: '+filename+'.geo\n')
        caseOut.close()
        meshOut=open(filename+'.geo','w')
        meshOut.write('Ensight Gold\n')
        meshOut.write('Unstructured Tetrahedral Mesh\n')
        meshOut.write('node id given\n')
        meshOut.write('element id given\n')
        #extents = 'extents\n %12.5E %12.5E\n %12.5E %12.5E\n %12.5E %12.5E\n' % (self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax)
        #meshOut.write('extents\n'+`self.xmin`+' '+`self.xmax`+'\n')
        meshOut.write('part \n'+'%10i\n' % 1)
        if description:
            meshOut.write(description+'\n')
        else:
            meshOut.write('A Mesh\n')
        meshOut.write('coordinates\n'+'%10i\n' % self.nNodes_global)
        for nN in range(self.nNodes_global):
            ensightNodeNumber = (nN+base)
            meshOut.write('%10i\n' % ensightNodeNumber)
        for nN in range(self.nNodes_global):
            meshOut.write('%12.5E\n' % self.nodeArray[nN,0])
        for nN in range(self.nNodes_global):
            meshOut.write('%12.5E\n' % self.nodeArray[nN,1])
        for nN in range(self.nNodes_global):
            meshOut.write('%12.5E\n' % self.nodeArray[nN,2])
        meshOut.write('tetra4\n'+'%10i\n' % self.nTetrahedra_global)
        for TN in range(self.nTetrahedra_global):
            ensightElementNumber = TN + base
            meshOut.write('%10i\n' % ensightElementNumber)
        TA = self.tetrahedronArray
        for TN in range(self.nTetrahedra_global):
            meshOut.write('%10i%10i%10i%10i\n' % (TA[TN,0]+base,
                                                  TA[TN,1]+base,
                                                  TA[TN,2]+base,
                                                  TA[TN,3]+base))
        meshOut.close()

    def writeBoundaryMeshEnsight(self,filename,description=None):
        base=1
        #write the casefile
        caseOut=open(filename+'Boundary.case','w')
        caseOut.write('FORMAT\n'+'type: ensight gold\n')
        caseOut.write('GEOMETRY\n'+'model: '+filename+'Boundary.geo\n')
        caseOut.close()
        meshOut=open(filename+'Boundary.geo','w')
        meshOut.write('Unstructured Triangular Surface Mesh\n\n')
        meshOut.write('node id given\n')
        meshOut.write('element id given\n')
        meshOut.write('part \n'+'%10i\n' % 1)
        if description:
            meshOut.write(description+'\n')
        else:
            meshOut.write('A Mesh\n')
        meshOut.write('coordinates\n'+'%10i\n' % self.nExteriorNodes_global)
        for nN in range(self.nExteriorNodes_global):
            ensightNodeNumber = (nN+base)
            meshOut.write('%10i\n' % ensightNodeNumber)
        for nN in range(self.nExteriorNodes_global):
            meshOut.write('%12.5E\n' %
                          self.nodeArray[self.exteriorNodeArray[nN],0])
        for nN in range(self.nExteriorNodes_global):
            meshOut.write('%12.5E\n' %
                          self.nodeArray[self.exteriorNodeArray[nN],1])
        for nN in range(self.nExteriorNodes_global):
            meshOut.write('%12.5E\n' %
                          self.nodeArray[self.exteriorNodeArray[nN],2])
        meshOut.write('tria3\n'+'%10i\n' % self.nExteriorTriangles_global)
        for tN in range(self.nExteriorTriangles_global):
            ensightElementNumber = tN + base
            meshOut.write('%10i\n' % ensightElementNumber)
        tA = self.triangleArray
        for tN in self.exteriorTriangleArray:
            meshOut.write('%10i%10i%10i\n' %
                          (self.globalToExteriorNodeArray[tA[tN,0]]+base,
                           self.globalToExteriorNodeArray[tA[tN,1]]+base,
                           self.globalToExteriorNodeArray[tA[tN,2]]+base))
        meshOut.close()

    def writeMeshXdmf(self,ar,name='',t=0.0,init=False,meshChanged=False,Xdmf_ElementTopology="Tetrahedron",tCount=0):
        if self.arGridCollection is not None:
            init = False
        elif not init:
            self.arGridCollection = ar.domain.find("Grid")
        if init:
            self.arGridCollection = SubElement(ar.domain,"Grid",{"Name":"Mesh "+name,
                                                               "GridType":"Collection",
                                                               "CollectionType":"Temporal"})
        if self.arGrid is None or self.arTime.get('Value') != str(t):
            if ar.global_sync:
                #
                #topology and geometry
                #
                self.arGrid = SubElement(self.arGridCollection,"Grid",{"GridType":"Uniform"})
                self.arTime = SubElement(self.arGrid,"Time",{"Value":str(t),"Name":str(tCount)})
                topology = SubElement(self.arGrid,"Topology",
                                      {"Type":Xdmf_ElementTopology,
                                       "NumberOfElements":str(self.globalMesh.nElements_global)})
                elements = SubElement(topology,"DataItem",
                                      {"Format":ar.dataItemFormat,
                                       "DataType":"Int",
                                       "Dimensions":"%i %i" % (self.globalMesh.nElements_owned,
                                                               self.nNodes_element)})
                geometry = SubElement(self.arGrid,"Geometry",{"Type":"XYZ"})
                nodes    = SubElement(geometry,"DataItem",
                                      {"Format":ar.dataItemFormat,
                                       "DataType":"Float",
                                       "Precision":"8",
                                       "Dimensions":"%i %i" % (self.globalMesh.nNodes_global,3)})
                #material types
                elementMaterialTypes = SubElement(self.arGrid,"Attribute",{"Name":"elementMaterialTypes",
                                                                           "AttributeType":"Scalar",
                                                                           "Center":"Cell"})
                elementMaterialTypesValues = SubElement(elementMaterialTypes,"DataItem",
                                                        {"Format":ar.dataItemFormat,
                                                         "DataType":"Int",
                                                         "Dimensions":"%i" % (self.globalMesh.nElements_owned,)})
                if ar.hdfFile is not None:
                    if ar.has_h5py:
                        elements.text = ar.hdfFilename+":/elements"+name+`tCount`
                        nodes.text = ar.hdfFilename+":/nodes"+name+`tCount`
                        elementMaterialTypesValues.text = ar.hdfFilename+":/"+"elementMaterialTypes"+"_t"+str(tCount)
                        if init or meshChanged:
                            ar.create_dataset_sync('elements'+name+`tCount`,
                                                   offsets = self.globalMesh.elementOffsets_subdomain_owned,
                                                   data = self.globalMesh.nodeNumbering_subdomain2global[self.elementNodesArray[:self.nElements_owned]])
                            ar.create_dataset_sync('nodes'+name+`tCount`,
                                                   offsets = self.globalMesh.nodeOffsets_subdomain_owned,
                                                   data = self.nodeArray[:self.nNodes_owned])
                            ar.create_dataset_sync("elementMaterialTypes"+"_t"+str(tCount),
                                                   offsets = self.globalMesh.elementOffsets_subdomain_owned,
                                                   data = self.elementMaterialTypes[:self.nElements_owned])
                    else:
                        assert False, "global_sync  not supported with pytables"
                else:
                    assert False, "global_sync not supported  with text heavy data"
            else:
                #
                #topology and geometry
                #
                self.arGrid = SubElement(self.arGridCollection,"Grid",{"GridType":"Uniform"})
                self.arTime = SubElement(self.arGrid,"Time",{"Value":str(t),"Name":str(tCount)})
                topology = SubElement(self.arGrid,"Topology",
                                      {"Type":Xdmf_ElementTopology,
                                       "NumberOfElements":str(self.nElements_owned)})
                elements = SubElement(topology,"DataItem",
                                      {"Format":ar.dataItemFormat,
                                       "DataType":"Int",
                                       "Dimensions":"%i %i" % (self.nElements_owned,self.nNodes_element)})
                geometry = SubElement(self.arGrid,"Geometry",{"Type":"XYZ"})
                nodes    = SubElement(geometry,"DataItem",
                                      {"Format":ar.dataItemFormat,
                                       "DataType":"Float",
                                       "Precision":"8",
                                       "Dimensions":"%i %i" % (self.nNodes_global,3)})
                #material types
                elementMaterialTypes = SubElement(self.arGrid,"Attribute",{"Name":"elementMaterialTypes",
                                                                           "AttributeType":"Scalar",
                                                                           "Center":"Cell"})
                elementMaterialTypesValues = SubElement(elementMaterialTypes,"DataItem",
                                                        {"Format":ar.dataItemFormat,
                                                         "DataType":"Int",
                                                         "Dimensions":"%i" % (self.nElements_owned,)})
                if ar.hdfFile is not None:
                    if ar.has_h5py:
                        elements.text = ar.hdfFilename+":/elements"+`ar.comm.rank()`+name+`tCount`
                        nodes.text = ar.hdfFilename+":/nodes"+`ar.comm.rank()`+name+`tCount`
                        elementMaterialTypesValues.text = ar.hdfFilename+":/"+"elementMaterialTypes"+"_p"+`ar.comm.rank()`+"_t"+str(tCount)
                        if init or meshChanged:
                            ar.create_dataset_async('elements'+`ar.comm.rank()`+name+`tCount`, data = self.elementNodesArray[:self.nElements_owned])
                            ar.create_dataset_async('nodes'+`ar.comm.rank()`+name+`tCount`, data = self.nodeArray)
                            ar.create_dataset_async("elementMaterialTypes"+"_p"+`ar.comm.rank()`+"_t"+str(tCount), data = self.elementMaterialTypes[:self.nElements_owned])
                    else:
                        elements.text = ar.hdfFilename+":/elements"+name+`tCount`
                        nodes.text = ar.hdfFilename+":/nodes"+name+`tCount`
                        elementMaterialTypesValues.text = ar.hdfFilename+":/"+"elementMaterialTypes"+str(tCount)
                        if init or meshChanged:
                            ar.hdfFile.createArray("/",'elements'+name+`tCount`,self.elementNodesArray[:self.nElements_owned])
                            ar.hdfFile.createArray("/",'nodes'+name+`tCount`,self.nodeArray)
                            ar.hdfFile.createArray("/","elementMaterialTypes"+str(tCount),self.elementMaterialTypes[:self.nElements_owned])
                else:
                    SubElement(elements,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/elements"+name+".txt"})
                    SubElement(nodes,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/nodes"+name+".txt"})
                    SubElement(elementMaterialTypesValues,"xi:include",{"parse":"text","href":"./"+ar.textDataDir+"/"+"elementMaterialTypes"+str(tCount)+".txt"})
                    if init or meshChanged:
                        np.savetxt(ar.textDataDir+"/elements"+name+".txt",self.elementNodesArray[:self.nElements_owned],fmt='%d')
                        np.savetxt(ar.textDataDir+"/nodes"+name+".txt",self.nodeArray)
                        np.savetxt(ar.textDataDir+"/"+"elementMaterialTypes"+str(tCount)+".txt",self.elementMaterialTypes[:self.nElements_owned])

class MultilevelTetrahedralMesh(MultilevelMesh):
    """A hierarchical multilevel mesh with tetrahedral cells"""
    def __init__(self,
                 nx, ny, nz,
                 x=0.0, y=0.0, z=0.0,
                 Lx=1.0, Ly=1.0, Lz=1.0,
                 refinementLevels=1,
                 skipInit=False,
                 nLayersOfOverlap=1,
                 parallelPartitioningType=MeshParallelPartitioningTypes.node):
        import cmeshTools
        import Comm
        MultilevelMesh.__init__(self)
        self.useC = True
        self.nLayersOfOverlap = nLayersOfOverlap; self.parallelPartitioningType = parallelPartitioningType
        logEvent("Generating tetrahedral mesh")
        if not skipInit:
            if self.useC:
                self.meshList.append(TetrahedralMesh())
                self.meshList[0].generateTetrahedralMeshFromRectangularGrid(nx,ny,nz,Lx,Ly,Lz)
                self.cmultilevelMesh = cmeshTools.CMultilevelMesh(self.meshList[0].cmesh,refinementLevels)
                self.buildFromC(self.cmultilevelMesh)
                self.meshList[0].nodeArray[:,0] += x
                self.meshList[0].nodeArray[:,1] += y
                self.meshList[0].nodeArray[:,2] += z
                self.meshList[0].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)
                for l in range(1,refinementLevels):
                    self.meshList.append(TetrahedralMesh())
                    self.meshList[l].cmesh = self.cmeshList[l]
                    self.meshList[l].buildFromC(self.cmeshList[l])
                    self.meshList[l].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)
            else:
                grid=RectangularGrid(nx,ny,nz,Lx,Ly,Lz)
                self.meshList.append(TetrahedralMesh())
                self.meshList[0].rectangularToTetrahedral(grid)
                self.meshList[0].nodeArray[:,0] += x
                self.meshList[0].nodeArray[:,1] += y
                self.meshList[0].nodeArray[:,2] += z
                self.elementChildren=[]
                logEvent(self.meshList[0].meshInfo())
                for l in range(1,refinementLevels):
                    self.refine()
                    logEvent(self.meshList[-1].meshInfo())
                self.buildArrayLists()
    def generateFromExistingCoarseMesh(self,mesh0,refinementLevels,nLayersOfOverlap=1,
                                       parallelPartitioningType=MeshParallelPartitioningTypes.node):
        import cmeshTools
        #blow away or just trust garbage collection
        self.nLayersOfOverlap=nLayersOfOverlap;self.parallelPartitioningType=parallelPartitioningType
        self.meshList = []
        self.elementParents = None
        self.cmultilevelMesh = None
        if self.useC:
            self.meshList.append(mesh0)
            logEvent("cmeshTools.CMultilevelMesh")
            self.cmultilevelMesh = cmeshTools.CMultilevelMesh(self.meshList[0].cmesh,refinementLevels)
            logEvent("buildFromC")
            self.buildFromC(self.cmultilevelMesh)
            logEvent("partitionMesh")
            self.meshList[0].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)
            for l in range(1,refinementLevels):
                self.meshList.append(TetrahedralMesh())
                self.meshList[l].cmesh = self.cmeshList[l]
                self.meshList[l].buildFromC(self.meshList[l].cmesh)
                self.meshList[l].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)
        else:
            grid=RectangularGrid(nx,ny,nz,Lx,Ly,Lz)
            self.meshList.append(TetrahedralMesh())
            self.meshList[0].rectangularToTetrahedral(grid)
            self.meshList[0].subdomainMesh = self.meshList[0]
            self.elementChildren=[]
            logEvent(self.meshList[0].meshInfo())
            for l in range(1,refinementLevels):
                self.refine()
                self.meshList[l].subdomainMesh = self.meshList[l]
                logEvent(self.meshList[-1].meshInfo())
            self.buildArrayLists()

    def generatePartitionedMeshFromPUMI(self,mesh0,refinementLevels,nLayersOfOverlap=1):
        import cmeshTools
        self.meshList = []
        self.meshList.append(mesh0)
        self.cmultilevelMesh = cmeshTools.CMultilevelMesh(self.meshList[0].cmesh,refinementLevels)
        self.buildFromC(self.cmultilevelMesh)
        self.elementParents = None
        self.elementChildren=[]

    def generatePartitionedMeshFromTetgenFiles(self,filebase,base,mesh0,refinementLevels,nLayersOfOverlap=1,
                                               parallelPartitioningType=MeshParallelPartitioningTypes.node):
        import cmeshTools
        assert(refinementLevels==1)
        assert(parallelPartitioningType==MeshParallelPartitioningTypes.node)
        assert(nLayersOfOverlap<=1)
        mesh0.cmesh = cmeshTools.CMesh()
        #blow away or just trust garbage collection
        self.nLayersOfOverlap=nLayersOfOverlap;self.parallelPartitioningType=parallelPartitioningType
        self.meshList = []
        self.elementParents = None
        self.cmultilevelMesh = None
        self.meshList.append(mesh0)
        logEvent("cmeshTools.CMultilevelMesh")
        self.cmultilevelMesh = cmeshTools.CMultilevelMesh(self.meshList[0].cmesh,refinementLevels)
        logEvent("buildFromC")
        self.buildFromC(self.cmultilevelMesh)
        logEvent("partitionMesh")
        self.meshList[0].partitionMeshFromFiles(filebase,base,nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)
    def refine(self):
        self.meshList.append(TetrahedralMesh())
        childrenDict = self.meshList[-1].refine(self.meshList[-2])
        self.elementChildren.append(childrenDict)
    def computeGeometricInfo(self):
        for m in self.meshList:
            m.computeGeometricInfo()

class MultilevelHexahedralMesh(MultilevelMesh):
    """A hierarchical multilevel mesh with hexahedral cells"""
    def __init__(self,
                 nx, ny, nz,
                 px=0, py=0, pz=0,
                 x=0.0, y=0.0, z=0.0,
                 Lx=1.0, Ly=1.0, Lz=1.0,
                 refinementLevels=1,
                 skipInit=False,
                 nLayersOfOverlap=1,
                 parallelPartitioningType=MeshParallelPartitioningTypes.node):
        import cmeshTools
        import Comm
        MultilevelMesh.__init__(self)
        if refinementLevels == 1:
            self.useC = True
        else:
            self.useC = False
        self.nLayersOfOverlap = nLayersOfOverlap; self.parallelPartitioningType = parallelPartitioningType
        logEvent("Generating hexahedral mesh")
        if not skipInit:
            if self.useC:
                self.meshList.append(HexahedralMesh())
                self.meshList[0].generateHexahedralMeshFromRectangularGrid(nx,ny,nz,Lx,Ly,Lz)
                self.meshList[0].nodeArray[:,0] += x
                self.meshList[0].nodeArray[:,1] += y
                self.meshList[0].nodeArray[:,2] += z
                self.cmultilevelMesh = cmeshTools.CMultilevelMesh(self.meshList[0].cmesh,refinementLevels)
                self.buildFromC(self.cmultilevelMesh)
                self.meshList[0].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)
                for l in range(1,refinementLevels):
                    self.meshList.append(HexahedralMesh())
                    self.meshList[l].cmesh = self.cmeshList[l]
                    self.meshList[l].buildFromC(self.cmeshList[l])
                    self.meshList[l].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)
            else:
                grid=RectangularGrid(nx,ny,nz,Lx,Ly,Lz)
                self.meshList.append(HexahedralMesh())
                self.elementChildren=[]
                self.meshList[0].sigmaMax=0.0
                logEvent(self.meshList[0].meshInfo())
                for l in range(1,refinementLevels):
                    self.refine()
                    self.meshList[-1].sigmaMax=0.0
                    logEvent(self.meshList[-1].meshInfo())
                self.buildArrayLists()
    def generateFromExistingCoarseMesh(self,mesh0,refinementLevels,nLayersOfOverlap=1,
                                       parallelPartitioningType=MeshParallelPartitioningTypes.node):
        import cmeshTools
        #blow away or just trust garbage collection
        self.nLayersOfOverlap=nLayersOfOverlap;self.parallelPartitioningType=parallelPartitioningType
        self.meshList = []
        self.elementParents = None
        self.cmultilevelMesh = None

        self.meshList.append(mesh0)
        self.cmultilevelMesh = cmeshTools.CMultilevelMesh(self.meshList[0].cmesh,refinementLevels)
        self.buildFromC(self.cmultilevelMesh)
        self.meshList[0].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)
        for l in range(1,refinementLevels):
            self.meshList.append(HexahedralMesh())
            self.meshList[l].cmesh = self.cmeshList[l]
            self.meshList[l].buildFromC(self.meshList[l].cmesh)
            self.meshList[l].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)


    def refine(self):
        raise NotImplementedError
        self.meshList.append(HexahedralMesh())
        childrenDict = self.meshList[-1].refine(self.meshList[-2])
        self.elementChildren.append(childrenDict)
    def computeGeometricInfo(self):
        for m in self.meshList:
            m.computeGeometricInfo()

def buildReferenceSimplex(nd=2):
    """
    Create and return a Proteus mesh object for the reference 
    element.

    Parameters
    ----------
    nd : int
        Dimension of reference element

    Returns
    -------
    mesh : :class:`proteus.MeshTools.TriangularMesh`
        Simplex mesh
    """
    from proteus import Domain
    from proteus import TriangleTools

    assert(nd in [1,2,3])

    if nd==1:
        pass # Note sure what needs to go here?!
    
    unit_simplex_domain = Domain.unitSimplex(nd)
    polyfile = "reference_element"
    unit_simplex_domain.writePoly(polyfile)

    if nd==2:
        tmesh = TriangleTools.TriangleBaseMesh(baseFlags="Yp",
                                               nbase=1,
                                               verbose=False)
        tmesh.readFromPolyFile(polyfile)
        mesh = tmesh.convertToProteusMesh(verbose=0)
        mesh.partitionMesh()
        mesh.globalMesh = mesh
        return mesh
    if nd==3:
        runTetgen(polyfile,
                  "Yp")
        mesh = genMeshWithTetgen(polyfile,
                                 nbase = 1)
        return mesh

class TriangularMesh(Mesh):
    """A mesh of triangles

    The nodes, edges, and triangles are indexed by their
    node tuples. The corresponding lists are derived from the dictionaries, and
    sorted lexicographically. The global node numbers are redefined to
    give a lexicographic ordering.

    The mesh can be generated from a rectangular grid and refined using either
    3t or Freudenthal-Bey global refinement.
    """

    def __init__(self):
        Mesh.__init__(self)
        self.nodeDict={}
        self.edgeDict={}
        self.triangleDict={}
        self.triangleList=[]
        self.oldToNewNode=[]
    def meshType(self):
        return 'simplex'
    def computeGeometricInfo(self):
        import cmeshTools
        cmeshTools.computeGeometricInfo_triangle(self.cmesh)
    def generateTriangularMeshFromRectangularGrid(self,nx,ny,Lx,Ly,triangleFlag=1):
        import cmeshTools
        self.cmesh = cmeshTools.CMesh()
        cmeshTools.generateTriangularMeshFromRectangularGrid(nx,ny,Lx,Ly,self.cmesh,triangleFlag)
        cmeshTools.allocateGeometricInfo_triangle(self.cmesh)
        cmeshTools.computeGeometricInfo_triangle(self.cmesh)
        self.buildFromC(self.cmesh)
    def rectangularToTriangularOriented(self,grid):
        #copy the nodes from the rectangular mesh
        #I want to be able to renumber latter without
        #changing the grid nodes, so I do deep copies here
        self.nodeList = [Node(n.N,n.p[X],n.p[Y],n.p[Z]) for n in grid.nodeList]
        self.nodeDict = dict([(n,n) for n in self.nodeList])
        for i in range(grid.nHx):
            for j in range(grid.nHy):
                k=0
                n0 = self.nodeList[grid.getNodeNumber(i,j,k)]
                n1 = self.nodeList[grid.getNodeNumber(i,j+1,k)]
                n2 = self.nodeList[grid.getNodeNumber(i+1,j,k)]
                n3 = self.nodeList[grid.getNodeNumber(i+1,j+1,k)]
                self.newTriangle([n0,n1,n3])
                self.newTriangle([n0,n2,n3])
        self.finalize()
        #self.buildListsEdges()
        #self.buildListsTriangles()
    def rectangularToTriangularOrientedOtherWay(self,grid):
        #copy the nodes from the rectangular mesh
        #I want to be able to renumber latter without
        #changing the grid nodes, so I do deep copies here
        self.nodeList = [Node(n.N,n.p[X],n.p[Y],n.p[Z]) for n in grid.nodeList]
        self.nodeDict = dict([(n,n) for n in self.nodeList])
        for i in range(grid.nHx):
            for j in range(grid.nHy):
                k=0
                n0 = self.nodeList[grid.getNodeNumber(i,j,k)]
                n1 = self.nodeList[grid.getNodeNumber(i,j+1,k)]
                n2 = self.nodeList[grid.getNodeNumber(i+1,j,k)]
                n3 = self.nodeList[grid.getNodeNumber(i+1,j+1,k)]
                self.newTriangle([n0,n2,n1])
                self.newTriangle([n2,n3,n1])
        self.finalize()
        #self.buildListsEdges()
        #self.buildListsTriangles()

    def rectangularToTriangularRedBlack(self,grid):
        #copy the nodes from the rectangular mesh
        #I want to be able to renumber latter without
        #changing the grid nodes, so I do deep copies here
        self.nodeList = [Node(n.N,n.p[X],n.p[Y],n.p[Z]) for n in grid.nodeList]
        self.nodeDict = dict([(n,n) for n in self.nodeList])
        self.triangleDict={}
        for i in range(grid.nHx):
            for j in range(grid.nHy):
                k=0
                #associate the element (i,j,k) with the
                #left, front, bottom node
                #get the left,front,bottom,node and its color
                if (grid.getColor(i,j,k) == grid.black):
                    b0 = self.nodeList[grid.getNodeNumber(i,j,k)]
                    r0 = self.nodeList[grid.getNodeNumber(i+1,j,k)]
                    r1 = self.nodeList[grid.getNodeNumber(i,j+1,k)]
                    b1 = self.nodeList[grid.getNodeNumber(i+1,j+1,k)]
                else:
                    r0 = self.nodeList[grid.getNodeNumber(i,j,k)]
                    b0 = self.nodeList[grid.getNodeNumber(i+1,j,k)]
                    b1 = self.nodeList[grid.getNodeNumber(i,j+1,k)]
                    r1 = self.nodeList[grid.getNodeNumber(i+1,j+1,k)]
                self.newTriangle([b0,r0,r1])
                self.newTriangle([b1,r0,r1])
        self.finalize()
        #self.buildListsEdges()
        #self.buildListsTriangles()
    #mwf debug switch to redblac
    rectangularToTriangular = rectangularToTriangularOrientedOtherWay#rectangularToTriangularOriented
    def generateFromTriangleMesh(self,ctrirep,base):
        import cmeshTools
        self.cmesh = cmeshTools.CMesh()
        cmeshTools.generateFromTriangleMesh(self.cmesh,ctrirep,base)
        cmeshTools.allocateGeometricInfo_triangle(self.cmesh)
        cmeshTools.computeGeometricInfo_triangle(self.cmesh)
        self.buildFromC(self.cmesh)
    def generateFromTriangleFiles(self,filebase,base):
        import cmeshTools
        self.cmesh = cmeshTools.CMesh()
        cmeshTools.generateFromTriangleFiles(self.cmesh,filebase,base)
        cmeshTools.allocateGeometricInfo_triangle(self.cmesh)
        cmeshTools.computeGeometricInfo_triangle(self.cmesh)
        self.buildFromC(self.cmesh)
    def writeTriangleFiles(self,filebase,base):
        import cmeshTools
        cmeshTools.writeTriangleFiles(self.cmesh,filebase,base)
    def generateFrom2DMFile(self,filebase,base=1):
        import cmeshTools
        self.cmesh = cmeshTools.CMesh()
        cmeshTools.generateFrom2DMFile(self.cmesh,filebase,base)
        cmeshTools.allocateGeometricInfo_triangle(self.cmesh)
        cmeshTools.computeGeometricInfo_triangle(self.cmesh)
        self.buildFromC(self.cmesh)

    def constructTriangularMeshOnRectangle(self,Lx,Ly,nx,ny,writeMesh=0,
                                           meshFileBase='mesh2d'):
        """
        wrapper function for making a triangular mesh on the rectangle
        [0,Lx] x [0,Ly].

        """
        nz = 1
        Lz = 1.0
        grid2d = RectangularGrid(nx,ny,nz,Lx,Ly,Lz)
        #grid2d.writeEdgesGnuplot('grid2d')
        #grid2d.viewMeshGnuplotPipe('grid2d')

        self.rectangularToTriangular(grid2d)

        if writeMesh == 1:
            #print mesh in gnuplot format
            self.writeEdgesGnuplot(meshFileBase)
            #can view with
            #self.viewMeshGnuplotPipe(meshFileBase)
        elif writeMesh == 2:
            self.writeEdgesMatlab(meshFileBase)
            #view in matlab with meshFileBase.m
        #end else

        return self

    def buildFromSets(self,triangleSet,edgeSet,nodeSet):
        self.nodeList = list(nodeSet)
        self.nodeDict = dict([(n,n) for n in self.nodeList])
        self.edgeList = list(edgeSet)
        self.edgeDict = dict([(e.nodes,e) for e in self.edgeList])
        self.triangleList = list(triangleSet)
        self.triangleDict = dict([(t.nodes,t) for t in self.triangleList])
        self.elementList = self.triangleList
        self.elementBoundaryList = self.edgeList
    def fixLocalNumbering(self):
        for tN in range(len(self.triangleList)):
            self.triangleList[tN].computeGeometricInfo()
            if edet(self.triangleList[tN].linearMap) < 0:
                newNodes = list(self.triangleList[tN].nodes)
                newNodes[2] = self.triangleList[tN].nodes[1]
                newNodes[1] = self.triangleList[tN].nodes[2]
                self.triangleList[tN].nodes = newNodes
    def finalize(self):
        self.buildLists()
        #self.fixLocalNumbering()
        self.buildArraysFromLists()
        #todo: build boundary mesh

    def buildLists(self):
        self.buildListsNodes()
        self.buildListsEdges()
        self.buildListsTriangles()
        self.elementList = self.triangleList
        self.elementBoundaryList = self.edgeList

    def buildListsNodes(self):
        keyList = self.nodeDict.keys()
        keyList.sort()
        self.nodeList=[]
        self.oldToNewNode=range(len(self.nodeDict))
        for nN,k in enumerate(keyList):
            self.oldToNewNode[self.nodeDict[k].N]=nN
            self.nodeDict[k].N = nN
            self.nodeList.append(self.nodeDict[k])

    def buildListsEdges(self):
        keyList = self.edgeDict.keys()
        keyList.sort()
        self.edgeList=[]
        for eN,k in enumerate(keyList):
            self.edgeDict[k].N = eN
            self.edgeList.append(self.edgeDict[k])

    def buildListsTriangles(self):
        keyList = self.triangleDict.keys()
        keyList.sort()
        self.triangleList=[]
        for tN,k in enumerate(keyList):
            self.triangleDict[k].N = tN
            self.triangleList.append(self.triangleDict[k])
        self.polygonList = self.triangleList

    def newTriangle(self,nodes):
        t = Triangle(len(self.triangleDict),nodes)
        self.triangleDict[t.nodes] = t
        self.registerEdges(t)
        return t

    def registerEdges(self,t):
        for en,e in enumerate(t.edges):
            if self.edgeDict.has_key(e.nodes):
                t.edges[en]=self.edgeDict[e.nodes]
            else:
                eN=len(self.edgeDict)
                e.N=eN
                self.edgeDict[e.nodes]=e

    def registerNode(self,node):
        if self.nodeDict.has_key(node):
            node = self.nodeDict[node]
        else:
            node.N = len(self.nodeDict)
            self.nodeDict[node] = node
        return node
    def buildLevelSetMesh(self,value,nodalValues):
        levelSetMesh = EdgeMesh()
        self.levelSetNodeNumbers = set()
        for t in self.triangleList:
            nodes={}
            for e in t.edges:
                nl = e.nodes[0]
                vl = nodalValues[nl.N]
                nr = e.nodes[1]
                vr = nodalValues[nr.N]
                if ((vl >= value and value >= vr) or
                    (vl <= value and value <= vr)):
                    if vl == vr:
                        newNl = Node(len(levelSetMesh.nodeDict),
                                     nl.p[X],
                                     nl.p[Y],
                                     nl.p[Z])
                        newNl = levelSetMesh.registerNode(newNl)
                        newNr = Node(len(levelSetMesh.nodeDict),
                                     nr.p[X],
                                     nr.p[Y],
                                     nr.p[Z])
                        newNr = levelSetMesh.registerNode(newNr)
                        levelSetMesh.newEdge([newNl,newNr])
                        self.levelSetNodeNumbers.add(nl.N)
                        self.levelSetNodeNumbers.add(nr.N)
                    elif value == vl:
                        newNode = Node(len(levelSetMesh.nodeDict),
                                       nl.p[X],
                                       nl.p[Y],
                                       nl.p[Z])
                        nodes[newNode] = newNode
                        self.levelSetNodeNumbers.add(nl.N)
                    elif value == vr and len(nodes) < 2:
                        newNode = Node(len(levelSetMesh.nodeDict),
                                       nr.p[X],
                                       nr.p[Y],
                                       nr.p[Z])
                        nodes[newNode] = newNode
                        self.levelSetNodeNumbers.add(nr.N)
                    else:
                        wr = (value - vl) / (vr - vl)
                        wl = (value - vr) / (vl - vr)
                        newPoint = nl.p*wl + nr.p*wr
                        newNode = Node(len(levelSetMesh.nodeDict),
                                       newPoint[X],
                                       newPoint[Y],
                                       newPoint[Z])
                        nodes[newNode] = newNode
                        self.levelSetNodeNumbers.add(nl.N)
                        self.levelSetNodeNumbers.add(nr.N)
                elif vl < value:
                    self.levelSetNodeNumbers.add(nl.N)
                elif vr < value:
                    self.levelSetNodeNumbers.add(nr.N)
            if len(nodes) == 0:
                pass
            elif len(nodes) == 1:
                print "singleton"
            elif len(nodes) == 2:
                newNodes=[]
                for n in nodes.values():
                    newNodes.append(levelSetMesh.registerNode(n))
                levelSetMesh.newEdge(newNodes)
            else:
                print "unexpected case in buildLevelSetMesh"
                print t.N
                for e in t.edges:
                    print e.N
                    for n in e.nodes:
                        print n.N
                        print n.p
                print "level set triangle"
                for n in nodes.values():
                    print n.p
        if len(levelSetMesh.edgeDict) == 0:
            print "level set does not cross any edges"
            return None
        else:
            levelSetMesh.finalize()
        return levelSetMesh
    def refine3t(self,oldMesh):
        childrenDict={}
        for t in oldMesh.triangleList:
            #deep copy old nodes because we'll renumber
            tNodes = [Node(eN,n.p[X],n.p[Y],n.p[Z])
                      for eN,n in enumerate(t.nodes)]
            for lnN,n in enumerate(tNodes): tNodes[lnN]=self.registerNode(n)
            #add new node
            t.computeGeometricInfo()
            newNode = Node(len(self.nodeDict),
                           t.barycenter[X],
                           t.barycenter[Y],
                           t.barycenter[Z])
            newNode = self.registerNode(newNode)
            t1=self.newTriangle([tNodes[0],tNodes[1],newNode])
            t2=self.newTriangle([tNodes[1],tNodes[2],newNode])
            t3=self.newTriangle([tNodes[2],tNodes[0],newNode])
            childrenDict[t.N]=[t1,t2,t3]
        self.finalize()
        return childrenDict

    def refineFreudenthalBey(self,oldMesh):
        logEvent("Refining the mesh using Freudenthal-Bey refinement")
        childrenDict={}
        for t in oldMesh.triangleDict.values():
            #deep copy old nodes because we'll renumber
            tNodes = [Node(nN,n.p[X],n.p[Y],n.p[Z])
                      for nN,n in enumerate(t.nodes)]
            for lnN,n in enumerate(tNodes): tNodes[lnN]=self.registerNode(n)
            #add new nodes (midpoints of edges)
            #use local edge tuples as keys
            newNodes={}
            for et,en in t.edgeMap.iteritems():
                t.edges[en].computeGeometricInfo()
                p = t.edges[en].barycenter
                newNodes[et] = Node(en,p[X],p[Y],p[Z])

            #set the global node numbers
            for k,n in newNodes.iteritems(): newNodes[k]=self.registerNode(n)
            #add corner triangles
            t1=self.newTriangle([tNodes[0],
                              newNodes[(0,1)],
                              newNodes[(0,2)]])
            t2=self.newTriangle([tNodes[1],
                              newNodes[(0,1)],
                              newNodes[(1,2)]])
            t3=self.newTriangle([tNodes[2],
                              newNodes[(0,2)],
                              newNodes[(1,2)]])
            #add center triangle
            t4=self.newTriangle([newNodes[(0,1)],
                              newNodes[(1,2)],
                              newNodes[(0,2)]])
            childrenDict[t.N]=[t1,t2,t3,t4]
        self.finalize()
        return childrenDict
        #for debugging: print each tet
        #self.edgeList=[]
        #Tlist = self.tetrahedronDict.values()
        #for T in Tlist:
        #    self.edgeList = self.edgeList + T.edges

    def refine(self,oldMesh):
        return self.refineFreudenthalBey(oldMesh)

    def meshInfo(self):
        minfo = """Number of triangles  : %d
Number of edges : %d
Number of nodes : %d\n""" % (self.nElements_global,
                             self.nElementBoundaries_global,
                             self.nNodes_global)
        if self.subdomainMesh != self:
            sinfo = self.subdomainMesh.meshInfo()
            info = "*** Global ***\n" + minfo + "\n*** Local ***\n" + sinfo
            return info
        return minfo

    def readMeshADH(self,filename,adhBase=1,suffix='3dm'):
        meshIn = open(filename+'.'+suffix,'r')
        firstLine = meshIn.readline()
        firstWords = firstLine.split()
        print "Reading object=%s from file=%s" % (firstWords[0],filename)
        line = meshIn.readline()
        columns = line.split()
        triangles = []
        triangleEdges=set()
        logEvent("Reading "+`filename`+ \
                " and building node lists for triangles, and edges")
        #assume triangles are ordered by triangle number
        while (columns[0] == 'E3T'):
            nodeNumbers = [int(c) - adhBase for c in columns[2:5]]
            nodeNumbers.sort()
            triangles.append(array.array('i',nodeNumbers))
            triangleEdges.update([(nodeNumbers[0],nodeNumbers[1]),
                                  (nodeNumbers[0],nodeNumbers[2]),
                                  (nodeNumbers[1],nodeNumbers[2])])
            line = meshIn.readline()
            columns = line.split()
        print "Building node list and dict"
        #assume nodes are ordered by node number
        while (len(columns) == 5):
            newNode = Node(int(columns[1]) - adhBase,
                           float(columns[2]),
                           float(columns[3]),
                           float(columns[4]))
            self.nodeList.append(newNode)
            self.nodeDict[newNode]=newNode
            line = meshIn.readline()
            columns = line.split()
        print "Number of triangles :"+`len(triangles)`
        print "Number of edges     :"+`len(triangleEdges)`
        print "Number of nodes     :"+`len(self.nodeList)`
        print "Number of objects   :"+\
              `len(triangleEdges)+len(triangles)+len(self.nodeList)`
        print "Building edge list"
        self.edgeList =[Edge(edgeNumber=eN,nodes=[self.nodeList[nN[0]],
                                                  self.nodeList[nN[1]]])
                        for eN,nN in enumerate(triangleEdges)]
        print "Building edge dict"
        self.edgeDict = dict([(e.nodes,e) for e in self.edgeList])
        print "Building triangle list"
        self.triangleList =[Triangle(triangleNumber=tN,
                                     nodes=[self.nodeList[nN[0]],
                                            self.nodeList[nN[1]],
                                            self.nodeList[nN[2]]],
                                     edgeDict=self.edgeDict)
                            for tN,nN in enumerate(triangles)]
        print "Building triangle dict"
        self.triangleDict = dict([(t.nodes,t) for t in self.triangleList])
        self.elementList = self.triangleList
        self.elementBoundaryList = self.edgeList

    def writeMeshXdmf(self,ar,name='',t=0.0,init=False,meshChanged=False,tCount=0,EB=False):
        Mesh.writeMeshXdmf(self,ar,name,t,init,meshChanged,"Triangle",tCount,EB=EB)
    def writeMeshEnsight(self,filename,description=None):
        base=1
        #write the casefile
        caseOut=open(filename+'.case','w')
        caseOut.write('FORMAT\n'+'type: ensight gold\n')
        caseOut.write('GEOMETRY\n'+'model: '+filename+'.geo\n')
        caseOut.close()
        meshOut=open(filename+'.geo','w')
        meshOut.write('Ensight Gold\n')
        meshOut.write('Unstructured Triangular Mesh\n')
        meshOut.write('node id given\n')
        meshOut.write('element id given\n')
        meshOut.write('part \n'+'%10i\n' % 1)
        if description:
            meshOut.write(description+'\n')
        else:
            meshOut.write('A Mesh\n')
        meshOut.write('coordinates\n'+'%10i\n' % self.nNodes_global)
        for nN in range(self.nNodes_global):
            meshOut.write('%10i\n' % (nN+base))
        for nN in range(self.nNodes_global):
            meshOut.write('%12.5E\n' % self.nodeArray[nN,0])
        for nN in range(self.nNodes_global):
            meshOut.write('%12.5E\n' % self.nodeArray[nN,1])
        for nN in range(self.nNodes_global):
            meshOut.write('%12.5E\n' % self.nodeArray[nN,2])
        meshOut.write('tria3\n'+'%10i\n' % self.nElements_global)
        for eN in range(self.nElements_global):
            meshOut.write('%10i\n' % (eN+base))
        for eN in range(self.nElements_global):
            meshOut.write('%10i%10i%10i\n' % tuple((nN+base) for nN in self.elementNodesArray[eN,:]))
        meshOut.close()

    def appendMeshEnsight(self,filename,description=None):
        base=1
        #write the casefile
        caseOut=open(filename+'.case','a')
        caseOut.write('GEOMETRY\n'+'model: '+filename+'.geo\n')
        caseOut.close()
        meshOut=open(filename+'.geo','w')
        meshOut.write('Unstructured Triangular Mesh\n\n')
        meshOut.write('node id given\n')
        meshOut.write('element id given\n')
        meshOut.write('part \n'+'%10i\n' % 1)
        if description:
            meshOut.write(description+'\n')
        else:
            meshOut.write('The whole mesh\n')
        meshOut.write('coordinates\n'+'%10i\n' % len(self.nodeList))
        for n in self.nodeList:
            nN = n.N+base
            meshOut.write('%10i\n' % nN)
        for n in self.nodeList:
            meshOut.write('%12.5E\n' % n.p[X])
        for n in self.nodeList:
            meshOut.write('%12.5E\n' % n.p[Y])
        for n in self.nodeList:
            meshOut.write('%12.5E\n' % n.p[Z])
        meshOut.write('tria3\n'+'%10i\n' % len(self.elementList))
        for e in self.elementList:
            eN = e.N + base
            meshOut.write('%10i\n' % eN)
        for e in self.elementList:
            meshOut.write('%10i%10i%10i\n' % tuple(n.N+base for n in e.nodes))
        meshOut.close()

    def writeMeshADH(self,filename,adhBase=1):
        import cmeshTools
        cmeshTools.write2dmFiles(self.cmesh,filename,adhBase)
    def writeAsymptote(self,fileprefix,L,x,units="m"):
        """
        Write a representation of the triangular mesh in the Asymptote vector graphics language
        """
        unitsize=4.0/L[0]
        f = open(fileprefix+".asy",'w')
        fileString="""
unitsize(4.0 inches / %(Lx)f);
size(5 inches);
real Lx=%(Lx)f;
real Ly=%(Ly)f;
real offset=0.0125Lx;
real x=%(x)f;
real y=%(y)f;
string strx="$%(Lx)2.2f\mbox{%(units)s}$";
string stry="$%(Ly)2.2f\mbox{%(units)s}$";
draw(strx,(x,y-offset)--(x+Lx,y-offset),S,black,Bars,Arrows,PenMargins);
draw(stry,(x-offset,y)--(x-offset,y+Ly),W,black,Bars,Arrows,PenMargins);
import graph;
import palette;
pen[] regionPens = Rainbow(NColors=%(nRegionFlags)d);
pen[] boundaryPens = Rainbow(NColors=%(nBoundaryFlags)d);
""" % {'Lx':L[0],'Ly':L[1],'x':x[0],'y':x[1],'units':units,
       'nRegionFlags':(max(self.elementMaterialTypes) - min(self.elementMaterialTypes)),
       'nBoundaryFlags':(max(self.elementBoundaryMaterialTypes)-min(self.elementBoundaryMaterialTypes))}
        #now draw triangles
        for t,tFlag in zip(self.elementNodesArray,self.elementMaterialTypes):
            fileString+="fill((%f,%f)--(%f,%f)--(%f,%f)--cycle,regionPens[%d]);\n" % (self.nodeArray[t[0]][0],self.nodeArray[t[0]][1],
                                                                                     self.nodeArray[t[1]][0],self.nodeArray[t[1]][1],
                                                                                     self.nodeArray[t[2]][0],self.nodeArray[t[2]][1],
                                                                                     tFlag-min(self.elementMaterialTypes))
        for eb,ebFlag in zip(self.elementBoundaryNodesArray,self.elementBoundaryMaterialTypes):
            if True:#ebFlag > 0:
                fileString+="draw((%f,%f)--(%f,%f),boundaryPens[%d]+linewidth(0.01));\n" % (self.nodeArray[eb[0]][0],self.nodeArray[eb[0]][1],
                                                                                           self.nodeArray[eb[1]][0],self.nodeArray[eb[1]][1],
                                                                                           ebFlag-min(self.elementBoundaryMaterialTypes))
        f.write(fileString)
        f.close()

#     def buildMatlabMeshDataStructures(self,meshFileBase='meshMatlab',writeToFile=True):
#         """
#         build array data structures for matlab finite element mesh representation
#         and write to a file to view and play with in matlatb

#         in matlab can then print mesh with

#         pdemesh(p,e,t)

#         where

#           p is the vertex or point matrix
#           e is the edge matrix, and
#           t is the element matrix

#         points matrix is [2 x num vertices]
#           format :
#              row 1 = x coord,
#              row 2 = y coord for nodes in mesh

#         edge matrix is [7 x num edges]
#           format:
#              row 1 = start vertex number
#              row 2 = end vertex number
#              row 3 = start value in edge parameterization, should be 0
#              row 4 = end   value in edge parameterization, should be 1
#              row 5 = global edge id, base 1
#              row 6 = subdomain on left? always 1 for now
#              row 7 = subdomain on right? always 0 for now

#         element matrix is [4 x num elements]
#             row 1 = vertex 1 global number
#             row 2 = vertex 2 global number
#             row 3 = vertex 3 global number
#             row 4 = triangle subdomain number
#          where 1,2,3 is a local counter clockwise numbering of vertices in
#            triangle

#          """
#         matlabBase = 1
#         p = np.zeros((2,self.nNodes_global),'d')
#         e = np.zeros((7,self.nElementBoundaries_global),'d')
#         t = np.zeros((4,self.nElements_global),'d')

#         #load p,e,t and write file
#         if writeToFile:
#             mfile = open(meshFileBase+'.m','w')
#         else:
#             mfile = open('/dev/null','w')
#         #
#         if writeToFile:
#             mfile.write('p = [ ... \n')
#         for nN in range(self.nNodes_global):
#             p[0,nN]=self.nodeArray[nN,0]
#             p[1,nN]=self.nodeArray[nN,1]
#             if writeToFile:
#                 mfile.write('%g %g \n' % tuple(p[:,nN]))
#         if writeToFile:
#             mfile.write(']; \n')
#             mfile.write("p = p\';\n")  #need transpose for matlab

#         if writeToFile:
#             mfile.write('e = [ ... \n')
#         for ebN in range(self.nElementBoundaries_global):
#             e[0,ebN]=self.elementBoundaryNodesArray[ebN,0] + matlabBase #global node number of start node base 1
#             e[1,ebN]=self.elementBoundaryNodesArray[ebN,1] + matlabBase #global node number of end node base 1
#             e[2,ebN]=0.0 #edge param. is 0 to 1
#             e[3,ebN]=1.0
#             e[4,ebN]=ebN + matlabBase  #global edge number base 1
#             e[5,ebN]=0 #subdomain to left
#             e[6,ebN]=1 #subdomain to right
#             if writeToFile:
#                 mfile.write('%g %g %g %g %g %g %g \n' % tuple(e[:,ebN]))
#         if writeToFile:
#             mfile.write(']; \n')
#             mfile.write("e = e\';\n")  #need transpose for matlab

#         #write triangles last
#         if writeToFile:
#             mfile.write('t = [ ... \n')
#         for eN in range(self.nElements_global):
#             t[0,eN]=self.elementNodesArray[eN,0]+matlabBase    #global node number for vertex 0
#             t[1,eN]=self.elementNodesArray[eN,1]+matlabBase    #global node number for vertex 0
#             t[2,eN]=self.elementNodesArray[eN,2]+matlabBase    #global node number for vertex 0
#             t[3,eN]=1                     #subdomain id
#             if writeToFile:
#                 mfile.write('%g %g %g %g \n' % tuple(t[:,eN]))
#         if writeToFile:
#             mfile.write(']; \n');
#             mfile.write("t = t\';\n") #need transpose for matlab



class QuadrilateralMesh(Mesh):
    """A mesh of quads

    The nodes, edges, and triangles are indexed by their
    node tuples. The corresponding lists are derived from the dictionaries, and
    sorted lexicographically. The global node numbers are redefined to
    give a lexicographic ordering.

    The mesh can be generated from a rectangular grid and refined using either
    3t or Freudenthal-Bey global refinement.
    """

    def __init__(self):
        Mesh.__init__(self)
        self.nodeDict={}
        self.edgeDict={}
        self.quadDict={}
        self.quadList=[]
        self.oldToNewNode=[]
        # tempoaray
        self.max_nNodeNeighbors_node = 4


    def buildFromSets(self,faceSet,edgeSet,nodeSet):
        self.nodeList = list(nodeSet)
        self.nodeDict = dict([(n,n) for n in self.nodeList])
        self.edgeList = list(edgeSet)
        self.edgeDict = dict([(e.nodes,e) for e in self.edgeList])
        self.quadList = list(faceSet)
        self.quadDict = dict([(t.nodes,t) for t in self.faceList])
        self.elementList = self.triangleList
        self.elementBoundaryList = self.edgeList

    def rectangularToQuadrilateral(self,grid,x=0.0,y=0.0,z=0.0):
        ''' WIP - I think this is the first function that needs to be
            written so that MultilevelQuadrilateralMesh can work.  This
            function does not call C functions.
        '''
        self.nodeList = [Node(n.N,n.p[X]+x,n.p[Y]+y,n.p[Z]+z) for n in grid.nodeList]
        # Is the following line necessary?
        self.nodeDict = dict([(n,n) for n in self.nodeList])
        for i in range(grid.nHx):
            for j in range(grid.nHy):
                k=0
                n0 = self.nodeList[grid.getNodeNumber(i,j,k)]
                n1 = self.nodeList[grid.getNodeNumber(i,j+1,k)]
                n2 = self.nodeList[grid.getNodeNumber(i+1,j+1,k)]
                n3 = self.nodeList[grid.getNodeNumber(i+1,j,k)]
                e0 = Edge(nodes=[n0,n1])
                e1 = Edge(nodes=[n1,n2])
                e2 = Edge(nodes=[n2,n3])
                e3 = Edge(nodes=[n3,n0])
                self.newQuadrilateral([e0,e1,e2,e3])
        self.finalize()

        
    def generateFromQuadFileIFISS(self,meshfile):
        ''' WIP - read a matlab.mat file containing IFISS vertices
        and elements
        '''
        import scipy.io
        griddata = scipy.io.loadmat(meshfile+'.mat')
        self.nodeList = [Node(nN,n[0],n[1],0.0) for nN,n in enumerate(griddata['vertices'])]
        # Is the following line necessary?
        self.nodeDict = dict([(n,n) for n in self.nodeList])
        for q in griddata['quads']:
            n0,n3,n2,n1 = q # clockwise ordering needed
            e0 = Edge(nodes=[self.nodeList[n0],self.nodeList[n1]])
            e1 = Edge(nodes=[self.nodeList[n1],self.nodeList[n2]])
            e2 = Edge(nodes=[self.nodeList[n2],self.nodeList[n3]])
            e3 = Edge(nodes=[self.nodeList[n3],self.nodeList[n0]])
            self.newQuadrilateral([e0,e1,e2,e3])
        self.finalize()
        for F,nN in griddata['bdyflags']:
            self.nodeMaterialTypes[nN] = F
        for ebNE in range(self.nExteriorElementBoundaries_global):
            ebN = self.exteriorElementBoundariesArray[ebNE]
            n0,n1 = self.elementBoundaryNodesArray[ebN]
            self.elementBoundaryMaterialTypes[ebN]=max(self.nodeMaterialTypes[n0],
                                                       self.nodeMaterialTypes[n1])


    def meshType(self):
        return 'cuboid'

    def meshInfo(self):
        minfo = """Number of quadrilaterals  : %d
Number of edges : %d
Number of nodes : %d\n""" % (self.nElements_global,
                             self.nElementBoundaries_global,
                             self.nNodes_global)
        if self.subdomainMesh != self:
            sinfo = self.subdomainMesh.meshInfo()
            info = "*** Global ***\n" + minfo + "\n*** Local ***\n" + sinfo
            return info
        return minfo

    def newQuadrilateral(self,edges):
        q = Quadrilateral(len(self.quadDict),edges)
        self.quadDict[q.nodes] = q
        self.registerEdges(q)
        return q

    def registerEdges(self,q):
        '''check if an edge is in the mesh dictionary
           if it is, point to existing entry
           otherwise, create a new entry
        '''
        for en,e in enumerate(q.edges):
            if self.edgeDict.has_key(e.nodes):
                q.edges[en]=self.edgeDict[e.nodes]
            else:
                eN=len(self.edgeDict)
                e.N=eN
                self.edgeDict[e.nodes]=e

    def registerNode(self,node):
        ''' check if a node is in the mesh dictionary
            if it is, point to existing entry
            otherwise, create a new entry
        '''
        if self.nodeDict.has_key(node):
            node = self.nodeDict[node]
        else:
            node.N = len(self.nodeDict)
            self.nodeDict[node] = node
        return node

    def refine(self,oldMesh):
        logEvent("Refining Using Standard Quadrilateral Refinement")
        import pdb
#        pdb.set_trace()
        childrenDict={}
        for q in oldMesh.quadDict.values():
            qNodes = [Node(nN,n.p[X],n.p[Y],n.p[Z]) for nN,n in enumerate(q.nodes)]
            for lnN,n in enumerate(qNodes): qNodes[lnN] = self.registerNode(n)
            q.computeGeometricInfo()

            newNodeLeft = Node(len(self.nodeDict),q.xmin,q.ymid,q.zmid)
            newNodeLeft = self.registerNode(newNodeLeft)
            newNodeTop = Node(len(self.nodeDict),q.xmid,q.ymax,q.zmid)
            newNodeTop = self.registerNode(newNodeTop)
            newNodeRight = Node(len(self.nodeDict),q.xmax,q.ymid,q.zmid)
            newNodeRight = self.registerNode(newNodeRight)
            newNodeBottom = Node(len(self.nodeDict),q.xmid,q.ymin,q.zmid)
            newNodeBottom = self.registerNode(newNodeBottom)
            newNodeMid = Node(len(self.nodeDict),q.xmid,q.ymid,q.zmid)
            newNodeMid = self.registerNode(newNodeMid)

            e1 = Edge(nodes=[qNodes[0],newNodeLeft])
            e2 = Edge(nodes=[newNodeLeft,newNodeMid])
            e3 = Edge(nodes=[newNodeMid,newNodeBottom])
            e4 = Edge(nodes=[newNodeBottom,qNodes[0]])
            e5 = Edge(nodes=[newNodeLeft,qNodes[1]])
            e6 = Edge(nodes=[qNodes[1],newNodeTop])
            e7 = Edge(nodes=[newNodeTop,newNodeMid])
            e8 = Edge(nodes=[newNodeTop,qNodes[2]])
            e9 = Edge(nodes=[qNodes[2], newNodeRight])
            e10 = Edge(nodes=[newNodeRight, newNodeMid])
            e11 = Edge(nodes=[qNodes[3],newNodeBottom])
            e12 = Edge(nodes=[newNodeRight,qNodes[3]])

            q1 = self.newQuadrilateral([e1,e2,e3,e4])
            self.registerEdges(q1)
            q2 = self.newQuadrilateral([e5,e6,e7,e2])
            self.registerEdges(q2)
            q3 = self.newQuadrilateral([e3,e10,e12,e11])
            self.registerEdges(q3)
            q4 = self.newQuadrilateral([e7,e8,e9,e10])
            self.registerEdges(q4)

            childrenDict[q.N]=[q1,q2,q3,q4]
        self.finalize()
        return childrenDict


    def generateQuadrialteralMeshFromRectangularGrid(self,nx,ny,Lx,Ly,triangleFlag=1):
        ''' WIP - This function needs to be constructed to allow MultilevelQuadrilateralMesh
            to run using C.  Implementing this will require some work to the mesh.cpp module. '''
        pass


    def finalize(self):
        ''' WIP '''
        self.buildLists()
        self.buildArraysFromLists()

    def buildLists(self):
        ''' WIP '''
        self.buildListsNodes()
        self.buildListsEdges()
        self.buildListsQuadrilaterals()
        self.elementList = self.quadList
        self.elementBoundaryList = self.edgeList

    def buildListsNodes(self):
        keyList = self.nodeDict.keys()
        keyList.sort()
        self.nodeList=[]
        self.oldToNewNode=range(len(self.nodeDict))
        for nN,k in enumerate(keyList):
            self.oldToNewNode[self.nodeDict[k].N]=nN
            self.nodeDict[k].N = nN
            self.nodeList.append(self.nodeDict[k])

    def buildListsEdges(self):
        keyList = self.edgeDict.keys()
        keyList.sort()
        self.edgeList=[]
        for eN,k in enumerate(keyList):
            self.edgeDict[k].N = eN
            self.edgeList.append(self.edgeDict[k])

    def buildListsQuadrilaterals(self):
        keyList = self.quadDict.keys()
        keyList.sort()
        self.quadList = []
        for qN,q in enumerate(keyList):
            self.quadDict[q].N = qN
            self.quadList.append(self.quadDict[q])
        self.polygonList = self.quadList

    def writeMeshXdmf(self,ar,name='',t=0.0,init=False,meshChanged=False,tCount=0,EB=False):
        Mesh.writeMeshXdmf(self,ar,name,t,init,meshChanged,"Quadrilateral",tCount,EB=EB)


class MultilevelTriangularMesh(MultilevelMesh):
    """A hierarchical  multilevel mesh of triangular cells"""
    import cmeshTools
    def __init__(self,
                 nx, ny, nz,
                 x=0.0, y=0.0, z=0.0,
                 Lx=1.0, Ly=1.0, Lz=1.0,
                 refinementLevels=1,
                 skipInit=False,
                 nLayersOfOverlap=1,
                 parallelPartitioningType=MeshParallelPartitioningTypes.node,triangleFlag=0):
        import cmeshTools
        MultilevelMesh.__init__(self)
        self.useC = True
        self.nLayersOfOverlap=nLayersOfOverlap; self.parallelPartitioningType = parallelPartitioningType
        #self.useC = False
        if not skipInit:
            if self.useC:
                self.meshList.append(TriangularMesh())
                self.meshList[0].generateTriangularMeshFromRectangularGrid(nx,ny,Lx,Ly,triangleFlag=triangleFlag)
                self.cmultilevelMesh = cmeshTools.CMultilevelMesh(self.meshList[0].cmesh,refinementLevels)
                self.buildFromC(self.cmultilevelMesh)
                self.meshList[0].nodeArray[:,0] += x
                self.meshList[0].nodeArray[:,1] += y
                self.meshList[0].nodeArray[:,2] += z
                self.meshList[0].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)
                for l in range(1,refinementLevels):
                    self.meshList.append(TriangularMesh())
                    self.meshList[l].cmesh = self.cmeshList[l]
                    self.meshList[l].buildFromC(self.meshList[l].cmesh)
                    self.meshList[l].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)
            else:
                grid=RectangularGrid(nx,ny,nz,Lx,Ly,Lz)
                self.meshList.append(TriangularMesh())
                self.meshList[0].rectangularToTriangular(grid)
                self.meshList[0].nodeArray[:,0] += x
                self.meshList[0].nodeArray[:,1] += y
                self.meshList[0].nodeArray[:,2] += z
                self.meshList[0].subdomainMesh = self.meshList[0]
                self.elementChildren=[]
                logEvent(self.meshList[0].meshInfo())
                for l in range(1,refinementLevels):
                    self.refine()
                    self.meshList[l].subdomainMesh = self.meshList[l]
                    logEvent(self.meshList[-1].meshInfo())
                self.buildArrayLists()
    #
    #mwf what's the best way to build from an existing mesh
    def generateFromExistingCoarseMesh(self,mesh0,refinementLevels,nLayersOfOverlap=1,
                                       parallelPartitioningType=MeshParallelPartitioningTypes.node):
        import cmeshTools
        #blow away or just trust garbage collection
        self.nLayersOfOverlap = nLayersOfOverlap; self.parallelPartitioningType = parallelPartitioningType
        self.meshList = []
        self.elementParents = None
        self.cmultilevelMesh = None
        if self.useC:
            self.meshList.append(mesh0)
            self.cmultilevelMesh = cmeshTools.CMultilevelMesh(self.meshList[0].cmesh,refinementLevels)
            self.buildFromC(self.cmultilevelMesh)
            self.meshList[0].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)
            for l in range(1,refinementLevels):
                self.meshList.append(TriangularMesh())
                self.meshList[l].cmesh = self.cmeshList[l]
                self.meshList[l].buildFromC(self.meshList[l].cmesh)
                self.meshList[l].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)
        else:
            grid=RectangularGrid(nx,ny,nz,Lx,Ly,Lz)
            self.meshList.append(TriangularMesh())
            self.meshList[0].rectangularToTriangular(grid)
            self.meshList[0].subdomainMesh = self.meshList[0]
            self.elementChildren=[]
            logEvent(self.meshList[0].meshInfo())
            for l in range(1,refinementLevels):
                self.refine()
                self.meshList[l].subdomainMesh = self.meshList[l]
                logEvent(self.meshList[-1].meshInfo())
            self.buildArrayLists()
    def generatePartitionedMeshFromPUMI(self,mesh0,refinementLevels,nLayersOfOverlap=1):
        import cmeshTools
        self.meshList = []
        self.meshList.append(mesh0)
        self.cmultilevelMesh = cmeshTools.CMultilevelMesh(self.meshList[0].cmesh,refinementLevels)
        self.buildFromC(self.cmultilevelMesh)
        self.elementParents = None
        self.elementChildren=[]

    def refine(self):
        self.meshList.append(TriangularMesh())
        childrenDict = self.meshList[-1].refine(self.meshList[-2])
        self.elementChildren.append(childrenDict)
        
    def computeGeometricInfo(self):
        for m in self.meshList:
            m.computeGeometricInfo()
    def locallyRefine(self,elementTagArray,flagForRefineType=0):
        """
        simple local refinement assuming elementTagArray[eN]=1 --> bisect

        flagForRefineType = 0 -- newest node, 1 -- 4T, 2 -- U4T
        """
        logEvent("MultilevelTriangularMesh:locallyRefine")
        if flagForRefineType == 0:
            logEvent("MultilevelTriangularMesh: calling cmeshTools.setNewestNodeBases")
            self.cmeshTools.setNewestNodeBases(2,self.cmultilevelMesh)
        if self.useC:
            logEvent("MultilevelTriangularMesh: calling locallRefineMultilevelMesh")
            self.cmeshTools.locallyRefineMultilevelMesh(2,self.cmultilevelMesh,elementTagArray,flagForRefineType)
            logEvent("MultilevelTriangularMesh: calling buildFromC")
            self.buildFromC(self.cmultilevelMesh)
            self.meshList.append(TriangularMesh())
            self.meshList[self.nLevels-1].cmesh = self.cmeshList[self.nLevels-1]
            self.meshList[self.nLevels-1].buildFromC(self.meshList[self.nLevels-1].cmesh)
            self.meshList[self.nLevels-1].partitionMesh(nLayersOfOverlap=self.nLayersOfOverlap,parallelPartitioningType=self.parallelPartitioningType)
        else:
            print """locallyRefine not implemented for self.useC= %s """ % (self.useC)
        #

class MultilevelQuadrilateralMesh(MultilevelMesh):
    """ A heirarchical multilevel mesh of quadrilaterals
       WIP """
    def __init__(self,
                 nx,ny,nz,
                 x=0.0,y=0.0,z=0.0,
                 Lx=1.0,Ly=1.0,Lz=1.0,
                 refinementLevels=1,
                 skipInit=False,
                 nLayersOfOverlap=1,
                 parallelPartitioningType=MeshParallelPartitioningTypes.node,triangleFlag=0):
        import cmeshTools
        MultilevelMesh.__init__(self)
        self.useC = False   # Implementing with C will take a bit more work. Disabling for now.
        self.nLayersOfOverlap=nLayersOfOverlap ; self.parallelPartitioningType = parallelPartitioningType
        if not skipInit:
            if self.useC:
                raise NotImplementedError ("C functionality sill not enabled for 2D quads")
            else:
                grid=RectangularGrid(nx,ny,nz,Lx,Ly,Lz)
                self.meshList.append(QuadrilateralMesh())
                self.meshList[0].rectangularToQuadrilateral(grid,x,y,z)
                self.meshList[0].subdomainMesh = self.meshList[0]
                self.elementChildren=[]
                logEvent(self.meshList[0].meshInfo())
                self.meshList[0].globalMesh = self.meshList[0]

                # The following four lines should be called elsewhere...Most of this is don in
                # the c-function calls that are not implemented yet for 2D quads
                self.meshList[0].nElements_owned = self.meshList[0].nElements_global
                self.meshList[0].nodeNumbering_subdomain2global.resize(self.meshList[0].nNodes_global)
                self.meshList[0].elementNumbering_subdomain2global.resize(self.meshList[0].nElements_global)
                self.meshList[0].nodeOffsets_subdomain_owned[-1] = self.meshList[0].nNodes_global
                self.meshList[0].nNodes_owned = self.meshList[0].nNodes_global
                self.meshList[0].elementOffsets_subdomain_owned[-1] = self.meshList[0].nElements_global

                for node in range(self.meshList[0].nNodes_global):
                    self.meshList[0].nodeNumbering_subdomain2global.itemset(node,node)
                for element in range(self.meshList[0].nElements_global):
                    self.meshList[0].elementNumbering_subdomain2global.itemset(element,element)
                for l in range(1,refinementLevels):
                    self.refine()
                    self.meshList[l].subdomainMesh = self.meshList[l]
                    logEvent(self.meshList[-1].meshInfo())
                self.buildArrayLists()

    def refine(self):
        self.meshList.append(QuadrilateralMesh())
        self.meshList[-1].globalMesh = self.meshList[-1]
        childrenDict = self.meshList[-1].refine(self.meshList[-2])

        # The following four lines should be called elsewhere...Most of this is don in
        # the c-function calls that are not implemented yet for 2D quads
        self.meshList[-1].nElements_owned = self.meshList[-1].nElements_global
        self.meshList[-1].nodeNumbering_subdomain2global.resize(self.meshList[-1].nNodes_global)
        self.meshList[-1].elementNumbering_subdomain2global.resize(self.meshList[-1].nElements_global)
        self.meshList[-1].nodeOffsets_subdomain_owned[-1] = self.meshList[-1].nNodes_global
        self.meshList[-1].nNodes_owned = self.meshList[-1].nNodes_global
        self.meshList[-1].elementOffsets_subdomain_owned[-1] = self.meshList[-1].nElements_global

        for node in range(self.meshList[-1].nNodes_global):
            self.meshList[-1].nodeNumbering_subdomain2global.itemset(node,node)
        for element in range(self.meshList[-1].nElements_global):
            self.meshList[-1].elementNumbering_subdomain2global.itemset(element,element)
        self.elementChildren.append(childrenDict)


class InterpolatedBathymetryMesh(MultilevelTriangularMesh):
    """A triangular mesh that interpolates bathymetry from a point cloud"""
    def __init__(self,
                 domain,
                 triangleOptions,
                 atol=1.0e-4,
                 rtol=1.0e-4,
                 maxElementDiameter=None,
                 maxLevels=20,
                 maxNodes=100000,
                 bathyType="points",#"grid"
                 bathyAssignmentScheme="interpolation",#"localAveraging","L2-projection","H1-projection"
                 errorNormType="L2", #L1,Linfty
                 refineType=0,
                 ):
        from matplotlib import tri as mpl_tri
        from scipy import interpolate as scipy_interpolate
        import TriangleTools
        if maxElementDiameter:
            self.maxElementDiameter = maxElementDiameter
        else:
            self.maxElementDiameter = np.inf
        self.atol = atol
        self.rtol = rtol
        self.maxLevels=maxLevels
        self.maxNodes=maxNodes
        self.domain = domain
        self.triangleOptions = triangleOptions
        self.bathyType=bathyType
        self.bathyAssignmentScheme=bathyAssignmentScheme
        self.errorNormType = errorNormType

        logEvent("InterpolatedBathymetryMesh: Calling Triangle to generate 2D coarse mesh for "+self.domain.name)
        tmesh = TriangleTools.TriangleBaseMesh(baseFlags=self.triangleOptions,
                                               nbase=1,
                                               verbose=10)
        tmesh.readFromPolyFile(domain.polyfile)

        logEvent("InterpolatedBathymetryMesh: Converting to Proteus Mesh")
        self.coarseMesh=tmesh.convertToProteusMesh(verbose=1)
        MultilevelTriangularMesh.__init__(self,0,0,0,skipInit=True,nLayersOfOverlap=0,
                                          parallelPartitioningType=MeshParallelPartitioningTypes.node)
        self.generateFromExistingCoarseMesh(self.coarseMesh,1,
                                            parallelPartitioningType=MeshParallelPartitioningTypes.node)
        self.computeGeometricInfo()
        print self.meshList[-1].volume
        #allocate some arrays based on the bathymetry data
        logEvent("InterpolatedBathymetryMesh:Allocating data structures for bathymetry interpolation algorithm")
        if bathyType == "points":
            self.nPoints_global = self.domain.bathy.shape[0]
            self.pointElementsArray_old = -np.ones((self.nPoints_global,),'i')
            self.pointElementsArray = -np.ones((self.nPoints_global,),'i')
            self.pointNodeWeightsArray = np.zeros((self.nPoints_global,3),'d')
            self.bathyInterpolant = scipy_interpolate.LinearNDInterpolator(self.domain.bathy[:,:2],self.domain.bathy[:,2])
            self.bathyNearestNeighbor = scipy_interpolate.NearestNDInterpolator(self.domain.bathy[:,:2], self.domain.bathy[:,2])
        elif bathyType == "grid":
            self.nPoints_global = self.domain.bathy.shape[0]
            self.pointElementsArray_old = -np.ones((self.nPoints_global,),'i')
            self.pointElementsArray = -np.ones((self.nPoints_global,),'i')
            self.pointNodeWeightsArray = np.zeros((self.nPoints_global,3),'d')
            x = self.domain.bathy[:self.domain.bathyGridDim[1],0]
            y = self.domain.bathy[:self.domain.bathyGridDim[0]*self.domain.bathyGridDim[1]:self.domain.bathyGridDim[1],1]
            z = self.domain.bathy[:,2].reshape(self.domain.bathyGridDim).transpose()
            self.bathyInterpolant = scipy_interpolate.RectBivariateSpline(x,y,z,kx=1,ky=1)
            #self.bathyInterpolant = scipy_interpolate.interp2d(x,y,z)
        #
        logEvent("InterpolatedBathymetryMesh: Locating points on initial mesh")
        self.locatePoints_initial(self.meshList[-1])
        logEvent("InterpolatedBathymetryMesh:setting mesh bathymetry from data")
        self.setMeshBathymetry(self.meshList[-1])
        logEvent("InterpolatedBathymetryMesh: tagging elements for refinement")
        self.tagElements(self.meshList[-1])
        levels = 0
        error = 1.0;
        while error >= 1.0 and self.meshList[-1].nNodes_global < self.maxNodes and levels < self.maxLevels:
            levels += 1
            logEvent("InterpolatedBathymetryMesh: Locally refining, level = %i" % (levels,))
            self.locallyRefine(self.meshList[-1].elementTags,flagForRefineType=refineType)
            logEvent("InterpolatedBathymetryMesh: interpolating bathymetry from parent mesh to refined mesh")
            self.interpolateBathymetry()
            logEvent("InterpolatedBathymetryMesh: Locating points on child mesh")
            self.locatePoints_refined(self.meshList[-1])
            logEvent("InterpolatedBathymetryMesh: setting mesh bathmetry from data")
            self.setMeshBathymetry(self.meshList[-1])
            logEvent("InterpolatedBathymetryMesh: tagging elements for refinement")
            error = self.tagElements(self.meshList[-1])
            logEvent("InterpolatedBathymetryMesh: error = %f atol = %f rtol = %f number of elements tagged = %i" % (error,self.atol,self.rtol,self.meshList[-1].elementTags.sum()))

    def setMeshBathymetry(self,mesh):
        if self.bathyAssignmentScheme == "interpolation":
            self.setMeshBathymetry_interpolate(mesh)
        elif self.bathyAssignmentScheme == "localAveraging":
            self.setMeshBathymetry_localAveraging(mesh)
        elif self.bathyAssignmentScheme == "L2-projection":
            raise NotImplementedError
        elif self.bathyAssignmentScheme == "H1-projection":
            raise NotImplementedError

    def setMeshBathymetry_interpolate(self,mesh):
        if self.bathyType == 'grid':
            mesh.nodeArray[:,2] = self.bathyInterpolant.ev(mesh.nodeArray[:,0],mesh.nodeArray[:,1])
        else:
            mesh.nodeArray[:,2] = self.bathyInterpolant(mesh.nodeArray[:,0],mesh.nodeArray[:,1])
            nI = np.isnan(mesh.nodeArray[:,2])
            mesh.nodeArray[nI,2] = self.bathyNearestNeighbor(mesh.nodeArray[nI,0],mesh.nodeArray[nI,1])

    def setMeshBathymetry_localAveraging(self,mesh):
        """
        calculate the arithmetic mean bathymetry of points inside each triangle and then assign the area-weighted average of the element means to each node
        """
        from FemTools import AffineMaps,ReferenceSimplex,LinearOnSimplexWithNodalBasis
        interpolationSpace = LinearOnSimplexWithNodalBasis(nd=2)
        #maps = AffineMaps(mesh,interpolationSpace.referenceElement,interpolationSpace)
        #maps.useC = True
        #calculate mean element height for each element
        #uses arithmetic mean, so it assumes the "patch" associated with each point the same size (weight)
        mesh.elementMeanZ = np.zeros((mesh.nElements_global,),'d')
        for pN in range(self.nPoints_global):
            eN = self.pointElementsArray[pN]
            if eN >= 0:
                if mesh.nPoints_element[eN] > 0:
                    mesh.elementMeanZ[eN] += self.domain.bathy[pN,2]/float(mesh.nPoints_element[eN])
                    mesh.nodeArray[mesh.elementNodesArray[eN,0],2] = 0.0
                    mesh.nodeArray[mesh.elementNodesArray[eN,1],2] = 0.0
                    mesh.nodeArray[mesh.elementNodesArray[eN,2],2] = 0.0
        #now assign the mesh node bathmetry as an area weighted average of the element mean
        sumArray = mesh.nodeArray[:,2].copy()
        sumArray[:]=0.0
        for eN in range(mesh.nElements_global):
            if mesh.nPoints_element[eN] > 0:#only calculate a contribution if this element contains a point
                #calculate triangle area and assign weighted average of element means to node
                xiArray = np.zeros((2,),'d')
                #
                grad_psi = np.zeros((interpolationSpace.dim,
                                        interpolationSpace.referenceElement.dim),
                                       'd')
                dx = np.zeros((interpolationSpace.referenceElement.dim),
                                 'd')
                jacobian = np.zeros((interpolationSpace.referenceElement.dim,
                                        interpolationSpace.referenceElement.dim),
                                       'd')
                inverseJacobian = np.zeros((interpolationSpace.referenceElement.dim,
                                               interpolationSpace.referenceElement.dim),
                                              'd')
                for j in interpolationSpace.range_dim:
                    grad_psi[j,:] = interpolationSpace.basisGradients[j](xiArray)#evaluate at zero because we can (psi is linear)
                jacobian.flat[:]=0.0
                inverseJacobian.flat[:]=0.0
                for j in interpolationSpace.range_dim:
                    J = mesh.elementNodesArray[eN,j]
                    for m in interpolationSpace.referenceElement.range_dim:
                        for n in interpolationSpace.referenceElement.range_dim:
                            jacobian[m,n] += mesh.nodeArray[J,m]*grad_psi[j,n]
                J = mesh.elementNodesArray[eN,0]
                inverseJacobian = inv(jacobian)
                area = 0.5*det(jacobian)
                sumArray[mesh.elementNodesArray[eN,0]] += area/mesh.nodeSupportArray[mesh.elementNodesArray[eN,0]]
                sumArray[mesh.elementNodesArray[eN,1]] += area/mesh.nodeSupportArray[mesh.elementNodesArray[eN,1]]
                sumArray[mesh.elementNodesArray[eN,2]] += area/mesh.nodeSupportArray[mesh.elementNodesArray[eN,2]]
                mesh.nodeArray[mesh.elementNodesArray[eN,0],2] += area*mesh.elementMeanZ[eN]/mesh.nodeSupportArray[mesh.elementNodesArray[eN,0]]
                mesh.nodeArray[mesh.elementNodesArray[eN,1],2] += area*mesh.elementMeanZ[eN]/mesh.nodeSupportArray[mesh.elementNodesArray[eN,1]]
                mesh.nodeArray[mesh.elementNodesArray[eN,2],2] += area*mesh.elementMeanZ[eN]/mesh.nodeSupportArray[mesh.elementNodesArray[eN,2]]
        #cek debug
        #print "sum of a nodes element areas divided by node support shoudl be 1 ",sumArray
    def locatePoints(self,mesh):
        """
        locate the element containing each point

        this should only be used on very coarse meshes
        """
        from FemTools import AffineMaps,ReferenceSimplex,LinearOnSimplexWithNodalBasis
        interpolationSpace = LinearOnSimplexWithNodalBasis(nd=2)
        #maps = AffineMaps(mesh,interpolationSpace.referenceElement,interpolationSpace)
        #maps.useC = False
        #find the elements that contain bathymetry points and calculate:
        # - for each element, the number of bathmetry points in that element
        # - for each node, the total area of the nodes elements that containing bathmetry points
        # - the area of each element
        # - the total area covered by elements containing bathmetry points
        mesh.nPoints_element = np.zeros((mesh.nElements_global,),'i')
        mesh.nodeSupportArray = np.zeros((mesh.nNodes_global,),'d')
        mesh.area_element =  np.zeros((mesh.nElements_global,),'d')
        self.pointElementsArray[:] = -1
        self.totalArea = 0.0
        for eN in range(mesh.nElements_global):
            #map points to reference space and test if it lies in the reference triangle
            xiArray = np.zeros((2,),'d')
            xiArray[:] = 0.0
            #
            grad_psi = np.zeros((interpolationSpace.dim,
                                    interpolationSpace.referenceElement.dim),
                                   'd')
            dx = np.zeros((interpolationSpace.referenceElement.dim),
                             'd')
            jacobian = np.zeros((interpolationSpace.referenceElement.dim,
                                    interpolationSpace.referenceElement.dim),
                                   'd')
            inverseJacobian = np.zeros((interpolationSpace.referenceElement.dim,
                                           interpolationSpace.referenceElement.dim),
                                          'd')
            for j in interpolationSpace.range_dim:
                grad_psi[j,:] = interpolationSpace.basisGradients[j](xiArray[0])#evalute at zero because we can (psi is linear)
            jacobian.flat[:]=0.0
            inverseJacobian.flat[:]=0.0
            for j in interpolationSpace.range_dim:
                J = mesh.elementNodesArray[eN,j]
                for m in interpolationSpace.referenceElement.range_dim:
                    for n in interpolationSpace.referenceElement.range_dim:
                        jacobian[m,n] += mesh.nodeArray[J,m]*grad_psi[j,n]
            J = mesh.elementNodesArray[eN,0]
            inverseJacobian = inv(jacobian)
            area = 0.5*det(jacobian)
            mesh.area_element[eN] = area
            self.totalArea += area
            for pN in range(self.nPoints_global):#can optimize by skipping previously found points
                xiArray[:] = 0.0
                dx[:]=self.domain.bathy[pN,:2]
                for m in interpolationSpace.referenceElement.range_dim:
                    dx[m]-=mesh.nodeArray[J,m]
                for m in interpolationSpace.referenceElement.range_dim:
                    for n in interpolationSpace.referenceElement.range_dim:
                        xiArray[m] += inverseJacobian[m,n]*dx[n]
                #barycentric coordinates are non-negative so we're in this element
                if xiArray[0] >=0.0 and xiArray[1] >= 0.0 and 1.0 - xiArray[0] - xiArray[1] >= 0.0:
                    self.pointElementsArray[pN] = eN
                    self.pointNodeWeightsArray[pN,0] = interpolationSpace.basis[0](xiArray)
                    self.pointNodeWeightsArray[pN,1] = interpolationSpace.basis[1](xiArray)
                    self.pointNodeWeightsArray[pN,2] = interpolationSpace.basis[2](xiArray)
        #count the number of points inside each element
        for pN in range(self.nPoints_global):
            if self.pointElementsArray[pN] >= 0:
                mesh.nPoints_element[self.pointElementsArray[pN]] += 1
        #add up the support area for each node
        for eN in range(mesh.nElements_global):
            if mesh.nPoints_element[eN] > 0:
                mesh.nodeSupportArray[mesh.elementNodesArray[eN,0]] += mesh.area_element[eN]
                mesh.nodeSupportArray[mesh.elementNodesArray[eN,1]] += mesh.area_element[eN]
                mesh.nodeSupportArray[mesh.elementNodesArray[eN,2]] += mesh.area_element[eN]
    def locatePoints_refined(self,mesh):
        """
        locate the element containing each point

        this should only be used on very coarse meshes
        """
        from FemTools import AffineMaps,ReferenceSimplex,LinearOnSimplexWithNodalBasis
        interpolationSpace = LinearOnSimplexWithNodalBasis(nd=2)
        #maps = AffineMaps(mesh,interpolationSpace.referenceElement,interpolationSpace)
        #maps.useC = False
        #find the elements that contain bathymetry points and calculate:
        # - for each element, the number of bathmetry points in that element
        # - for each node, the total area of the nodes elements that containing bathmetry points
        # - the area of each element
        # - the total area covered by elements containing bathmetry points
        mesh.nPoints_element = np.zeros((mesh.nElements_global,),'i')
        mesh.nodeSupportArray = np.zeros((mesh.nNodes_global,),'d')
        mesh.area_element =  np.zeros((mesh.nElements_global,),'d')
        self.totalArea = 0.0
        self.pointElementsArray_old[:] = self.pointElementsArray
        self.pointElementsArray[:] = -1
        for pN in range(self.nPoints_global):
            eN_parent = self.pointElementsArray_old[pN]
            for eN in self.elementChildrenArrayList[-1][self.elementChildrenOffsetsList[-1][eN_parent]:self.elementChildrenOffsetsList[-1][eN_parent+1]]:
                xiArray = np.zeros((2,),'d')
                xiArray[:] = 0.0
                grad_psi = np.zeros((interpolationSpace.dim,
                                        interpolationSpace.referenceElement.dim),
                                       'd')
                dx = np.zeros((interpolationSpace.referenceElement.dim),
                                 'd')
                jacobian = np.zeros((interpolationSpace.referenceElement.dim,
                                        interpolationSpace.referenceElement.dim),
                                       'd')
                inverseJacobian = np.zeros((interpolationSpace.referenceElement.dim,
                                               interpolationSpace.referenceElement.dim),
                                              'd')
                for j in interpolationSpace.range_dim:
                    grad_psi[j,:] = interpolationSpace.basisGradients[j](xiArray[0])#evalute at zero because we can (psi is linear)
                jacobian.flat[:]=0.0
                inverseJacobian.flat[:]=0.0
                for j in interpolationSpace.range_dim:
                    J = mesh.elementNodesArray[eN,j]
                    for m in interpolationSpace.referenceElement.range_dim:
                        for n in interpolationSpace.referenceElement.range_dim:
                            jacobian[m,n] += mesh.nodeArray[J,m]*grad_psi[j,n]
                J = mesh.elementNodesArray[eN,0]
                inverseJacobian = inv(jacobian)
                area = 0.5*det(jacobian)
                mesh.area_element[eN] = area
                self.totalArea += area
                xiArray[:] = 0.0
                dx[:]=self.domain.bathy[pN,:2]
                for m in interpolationSpace.referenceElement.range_dim:
                    dx[m]-=mesh.nodeArray[J,m]
                for m in interpolationSpace.referenceElement.range_dim:
                    for n in interpolationSpace.referenceElement.range_dim:
                        xiArray[m] += inverseJacobian[m,n]*dx[n]
                #barycentric coordinates are non-negative so we're in this element
                if xiArray[0] >=0.0 and xiArray[1] >= 0.0 and 1.0 - xiArray[0] - xiArray[1] >= 0.0:
                    self.pointElementsArray[pN] = eN
                    self.pointNodeWeightsArray[pN,0] = interpolationSpace.basis[0](xiArray)
                    self.pointNodeWeightsArray[pN,1] = interpolationSpace.basis[1](xiArray)
                    self.pointNodeWeightsArray[pN,2] = interpolationSpace.basis[2](xiArray)
        #count the number of points inside each element
        for pN in range(self.nPoints_global):
            if self.pointElementsArray[pN] >= 0:
                mesh.nPoints_element[self.pointElementsArray[pN]] += 1
        #add up the support area for each node
        for eN in range(mesh.nElements_global):
            if mesh.nPoints_element[eN] > 0:
                mesh.nodeSupportArray[mesh.elementNodesArray[eN,0]] += mesh.area_element[eN]
                mesh.nodeSupportArray[mesh.elementNodesArray[eN,1]] += mesh.area_element[eN]
                mesh.nodeSupportArray[mesh.elementNodesArray[eN,2]] += mesh.area_element[eN]
    def locatePoints_initial(self,mesh):
        """
        locate the element containing each point

        first find the nearest node, then loop over that node's  elements
        """
        from scipy.spatial import cKDTree
        from FemTools import AffineMaps,ReferenceSimplex,LinearOnSimplexWithNodalBasis
        interpolationSpace = LinearOnSimplexWithNodalBasis(nd=2)
        #find the elements that contain bathymetry points and calculate:
        # - for each element, the number of bathmetry points in that element
        # - for each node, the total area of the nodes elements that contain bathmetry points
        # - the area of each element
        # - the total area covered by elements containing bathmetry points
        mesh.nPoints_element = np.zeros((mesh.nElements_global,),'i')
        mesh.nodeSupportArray = np.zeros((mesh.nNodes_global,),'d')
        mesh.area_element =  np.zeros((mesh.nElements_global,),'d')
        self.totalArea = 0.0
        self.pointElementsArray[:] = -1
        tree = cKDTree(mesh.nodeArray[:,:2])
        xiArray = np.zeros((2,),'d')
        grad_psi = np.zeros((interpolationSpace.dim,
                                interpolationSpace.referenceElement.dim),
                               'd')
        dx = np.zeros((interpolationSpace.referenceElement.dim),
                         'd')
        jacobian = np.zeros((interpolationSpace.referenceElement.dim,
                                interpolationSpace.referenceElement.dim),
                               'd')
        inverseJacobian = np.zeros((interpolationSpace.referenceElement.dim,
                                       interpolationSpace.referenceElement.dim),
                                      'd')
        for pN in range(self.nPoints_global):
            (distance,nN) = tree.query(self.domain.bathy[pN,:2])
            for eN in mesh.nodeElementsArray[mesh.nodeElementOffsets[nN]:mesh.nodeElementOffsets[nN+1]]:
                xiArray[:] = 0.0
                for j in interpolationSpace.range_dim:
                    grad_psi[j,:] = interpolationSpace.basisGradients[j](xiArray[0])#evalute at zero because we can (psi is linear)
                jacobian.flat[:]=0.0
                inverseJacobian.flat[:]=0.0
                for j in interpolationSpace.range_dim:
                    J = mesh.elementNodesArray[eN,j]
                    for m in interpolationSpace.referenceElement.range_dim:
                        for n in interpolationSpace.referenceElement.range_dim:
                            jacobian[m,n] += mesh.nodeArray[J,m]*grad_psi[j,n]
                J = mesh.elementNodesArray[eN,0]
                inverseJacobian = inv(jacobian)
                area = 0.5*det(jacobian)
                mesh.area_element[eN] = area
                xiArray[:] = 0.0
                dx[:]=self.domain.bathy[pN,:2]
                for m in interpolationSpace.referenceElement.range_dim:
                    dx[m]-=mesh.nodeArray[J,m]
                for m in interpolationSpace.referenceElement.range_dim:
                    for n in interpolationSpace.referenceElement.range_dim:
                        xiArray[m] += inverseJacobian[m,n]*dx[n]
                #if the barycentric coordinates are non-negative we're in this element
                if xiArray[0] >=0.0 and xiArray[1] >= 0.0 and 1.0 - xiArray[0] - xiArray[1] >= 0.0:
                    self.pointElementsArray[pN] = eN
                    self.pointNodeWeightsArray[pN,0] = interpolationSpace.basis[0](xiArray)
                    self.pointNodeWeightsArray[pN,1] = interpolationSpace.basis[1](xiArray)
                    self.pointNodeWeightsArray[pN,2] = interpolationSpace.basis[2](xiArray)
        self.totalArea += mesh.area_element.sum()
        #count the number of points inside each element
        for pN in range(self.nPoints_global):
            if self.pointElementsArray[pN] >= 0:
                mesh.nPoints_element[self.pointElementsArray[pN]] += 1
        #add up the support area for each node
        for eN in range(mesh.nElements_global):
            if mesh.nPoints_element[eN] > 0:
                mesh.nodeSupportArray[mesh.elementNodesArray[eN,0]] += mesh.area_element[eN]
                mesh.nodeSupportArray[mesh.elementNodesArray[eN,1]] += mesh.area_element[eN]
                mesh.nodeSupportArray[mesh.elementNodesArray[eN,2]] += mesh.area_element[eN]
    def interpolateBathymetry(self):
        """
        interpolate bathymetry for the refinement from the parent  mesh
        """
        from proteus.FemTools import C0_AffineLinearOnSimplexWithNodalBasis,DOFBoundaryConditions,MultilevelProjectionOperators
        mlMeshTemp = MultilevelMesh(levels=2)
        mlMeshTemp.meshList = self.meshList[-2:]
        mlMeshTemp.nLevels=2
        mlMeshTemp.cmeshList = self.cmeshList[-2:]
        mlMeshTemp.elementParentsArrayList = self.elementParentsArrayList[-2:]
        mlMeshTemp.elementChildrenArrayList = self.elementChildrenArrayList[-1:]
        mlMeshTemp.elementChildrenOffsetsList = self.elementChildrenOffsetsList[-1:]
        nd=2
        TrialSpaceTypeDict = {0:C0_AffineLinearOnSimplexWithNodalBasis}
        trialSpaceDictParent = dict([ (cj,TrialSpaceType(mlMeshTemp.meshList[0],nd)) for (cj,TrialSpaceType) in TrialSpaceTypeDict.iteritems()])
        trialSpaceDictChild = dict([ (cj,TrialSpaceType(mlMeshTemp.meshList[1],nd)) for (cj,TrialSpaceType) in TrialSpaceTypeDict.iteritems()])
        trialSpaceDictList  = [trialSpaceDictParent,trialSpaceDictChild]
        offsetListList=[[0],[0]]
        strideListList=[[1],[1]]
        def getDBC(x,flag):
            return None
        bcDictList=[dict([(0,DOFBoundaryConditions(trialSpaceDictParent[0],getPointwiseBoundaryConditions=getDBC,weakDirichletConditions=False))]),
                    dict([(0,DOFBoundaryConditions(trialSpaceDictChild[0],getPointwiseBoundaryConditions=getDBC,weakDirichletConditions=False))])]
        self.meshTransfers = MultilevelProjectionOperators(
            mlMeshTemp,
            trialSpaceDictList,
            offsetListList,
            strideListList,
            bcDictList)
        zParent = self.meshList[-2].nodeArray[:,2].copy()
        zChild = self.meshList[-1].nodeArray[:,2].copy()
        self.meshTransfers.prolongList[-1].matvec(zParent,zChild)
        self.meshList[-1].nodeArray[:,2] = zChild
    def tagElements(self,mesh):
        """
        loop over points and calculate whether the interpolation error is within the tolerance

        this should only be used on very coarse meshes
        """
        mesh.elementTags = np.zeros((mesh.nElements_global,),'i')
        mesh.errorAverage_element =  np.zeros((mesh.nElements_global,),'d')
        errorInfty = 0.0
        mesh.elementTags[mesh.elementDiametersArray > self.maxElementDiameter ] = 1
        print mesh.elementTags
        for pN in range(self.nPoints_global):
            eN = self.pointElementsArray[pN]
            if eN >= 0:
                #print "pN = ",pN,"eN = ",eN,"nodes ",mesh.elementNodesArray[eN,:]
                zInterp = self.pointNodeWeightsArray[pN,0]*mesh.nodeArray[mesh.elementNodesArray[eN,0],2] +  \
                          self.pointNodeWeightsArray[pN,1]*mesh.nodeArray[mesh.elementNodesArray[eN,1],2] +  \
                          self.pointNodeWeightsArray[pN,2]*mesh.nodeArray[mesh.elementNodesArray[eN,2],2]
                errorPointwise = fabs(zInterp - self.domain.bathy[pN,2]) / (fabs(self.domain.bathy[pN,2])*self.rtol + self.atol)
                #print "error ",errorPointwise
                errorInfty = max(errorPointwise,errorInfty)
                mesh.errorAverage_element[eN] += (errorPointwise/float(mesh.nPoints_element[eN]))
                #print "error average",mesh.errorAverage_element[eN]
                if errorPointwise >= 1.0:
                    mesh.elementTags[eN] = 1
        if self.errorNormType == "L1":
            mesh.elementTags[:] = 0
            errorL1 = 0.0
            for eN in range(mesh.nElements_global):
                errorL1 += mesh.errorAverage_element[eN]*mesh.area_element[eN]
                if mesh.errorAverage_element[eN] >= 1.0:
                    mesh.elementTags[eN] = 1
            errorL1 /= self.totalArea#normalize by domain error to make error have units of length
            return errorL1
        if self.errorNormType == "L2":
            mesh.elementTags[:] = 0
            errorL2 = 0.0
            for eN in range(mesh.nElements_global):
                errorL2 += (mesh.errorAverage_element[eN])**2 * mesh.area_element[eN]
                if mesh.errorAverage_element[eN] >= 1.0:
                    mesh.elementTags[eN] = 1
            errorL2 = sqrt(errorL2)/self.totalArea#normalize by domain error to make error have units of length
            return errorL2
        else:
            #print "finished"
            return errorInfty

class EdgeMesh(Mesh):
    """A mesh of edges

    The nodes, and edges are indexed by their node tuples. The
    corresponding lists are derived from the dictionaries, and sorted
    lexicographically. The global node numbers are redefined to give a
    lexicographic ordering.
    """

    def __init__(self):
        Mesh.__init__(self)
        self.nodeDict={}
        self.edgeDict={}
        self.oldToNewNode=[]
    def computeGeometricInfo(self):
        import cmeshTools
        cmeshTools.computeGeometricInfo_edge(self.cmesh)
    def generateEdgeMeshFromRectangularGrid(self,nx,Lx):
        import cmeshTools
        self.cmesh = cmeshTools.CMesh()
        cmeshTools.generateEdgeMeshFromRectangularGrid(nx,Lx,self.cmesh)
        cmeshTools.allocateGeometricInfo_edge(self.cmesh)
        cmeshTools.computeGeometricInfo_edge(self.cmesh)
        self.buildFromC(self.cmesh)
        #mwf debug
        #print "EdgeMesh rect->edge after build nodes=%s " % (self.nodeArray)
    def rectangularToEdge(self,grid):
        #copy the nodes from the rectangular mesh
        #I want to be able to renumber latter without
        #changing the grid nodes, so I do deep copies here
        self.nodeList = [Node(n.N,n.p[X],n.p[Y],n.p[Z]) for n in grid.nodeList]
        self.nodeDict = dict([(n,n) for n in self.nodeList])
        for e in grid.edgeList:
            self.newEdge([self.nodeDict[e.nodes[0]],self.nodeDict[e.nodes[1]]])
        self.finalize()
        #self.buildListsEdges()

    def finalize(self):
        self.buildLists()
        self.buildArraysFromLists()
        #todo: build boundary mesh

    def buildLists(self):
        self.buildListsNodes()
        self.buildListsEdges()
        self.elementList = self.edgeList
        self.elementBoundaryList = self.nodeList

    def buildListsNodes(self):
        keyList = self.nodeDict.keys()
        keyList.sort()
        self.nodeList=[]
        self.oldToNewNode=range(len(self.nodeDict))
        for nN,k in enumerate(keyList):
            self.oldToNewNode[self.nodeDict[k].N]=nN
            self.nodeDict[k].N = nN
            self.nodeList.append(self.nodeDict[k])

    def buildListsEdges(self):
        keyList = self.edgeDict.keys()
        keyList.sort()
        self.edgeList=[]
        for eN,k in enumerate(keyList):
            self.edgeDict[k].N = eN
            self.edgeList.append(self.edgeDict[k])

    def newEdge(self,nodes):
        e = Edge(len(self.edgeDict),nodes)
        self.edgeDict[e.nodes] = e
        return e

    def registerNode(self,node):
        if self.nodeDict.has_key(node):
            node = self.nodeDict[node]
        else:
            node.N = len(self.nodeDict)
            self.nodeDict[node] = node
        return node

    def refine2e(self,oldMesh):
        childrenDict={}
        for e in oldMesh.edgeList:
            #deep copy old nodes because we'll renumber
            eNodes = [Node(eN,n.p[X],n.p[Y],n.p[Z])
                      for eN,n in enumerate(e.nodes)]
            for lnN,n in enumerate(eNodes): eNodes[lnN]=self.registerNode(n)
            #add new node
            e.computeGeometricInfo()
            newNode = Node(len(self.nodeDict),
                           e.barycenter[X],
                           e.barycenter[Y],
                           e.barycenter[Z])
            newNode = self.registerNode(newNode)
            e1=self.newEdge([eNodes[0],newNode])
            e2=self.newEdge([newNode,eNodes[1]])
            childrenDict[e.N]=[e1,e2]
        self.finalize()
        return childrenDict

    def refine(self,oldMesh):
        return self.refine2e(oldMesh)

    def meshInfo(self):
        minfo = """Number of edges : %d
Number of nodes : %d\n""" % (self.nElements_global,self.nNodes_global)
        if self.subdomainMesh != self:
            sinfo = self.subdomainMesh.meshInfo()
            info = "*** Global ***\n" + minfo + "\n*** Local ***\n" + sinfo
            return info
        return minfo

    def writeMeshADH(self,filename):
        pass
    def writeMeshXdmf(self,ar,name='',t=0.0,init=False,meshChanged=False,tCount=0):
        Mesh.writeMeshXdmf(self,ar,name,t,init,meshChanged,"Polyline",tCount)
    def writeMeshEnsight(self,filename,description=None):
        base=1
        #write the casefile
        caseOut=open(filename+'.case','w')
        caseOut.write('FORMAT\n'+'type: ensight gold\n')
        caseOut.write('GEOMETRY\n'+'model: '+filename+'.geo\n')
        caseOut.close()
        meshOut=open(filename+'.geo','w')
        meshOut.write('Ensight Gold\n')
        meshOut.write('Unstructured Edge Mesh\n')
        meshOut.write('node id given\n')
        meshOut.write('element id given\n')
        meshOut.write('part \n'+'%10i\n' % 1)
        if description:
            meshOut.write(description+'\n')
        else:
            meshOut.write('A Mesh\n')
        meshOut.write('coordinates\n'+'%10i\n' % self.nNodes_global)
        for nN in range(self.nNodes_global):
            meshOut.write('%10i\n' % (nN+base))
        for nN in range(self.nNodes_global):
            meshOut.write('%12.5E\n' % self.nodeArray[nN,0])
        for nN in range(self.nNodes_global):
            meshOut.write('%12.5E\n' % self.nodeArray[nN,1])
        for nN in range(self.nNodes_global):
            meshOut.write('%12.5E\n' % self.nodeArray[nN,2])
        meshOut.write('bar2\n'+'%10i\n' % self.nElements_global)
        for eN in range(self.nElements_global):
            meshOut.write('%10i\n' % (eN+base))
        for eN in range(self.nElements_global):
            meshOut.write('%10i%10i\n' % tuple((nN+base) for nN in self.elementNodesArray[eN,:]))
        meshOut.close()

class MultilevelEdgeMesh(MultilevelMesh):
    """A hierarchical multilevel mesh of intervals (edges)"""
    import cmeshTools
    def __init__(self,
                 nx, ny, nz,
                 x=0.0, y=0.0, z=0.0,
                 Lx=1.0, Ly=1.0, Lz=1.0,
                 refinementLevels=1,
                 nLayersOfOverlap=1,
                 parallelPartitioningType=MeshParallelPartitioningTypes.node):
        import cmeshTools
        MultilevelMesh.__init__(self)
        self.useC=True
        self.nLayersOfOverlap=nLayersOfOverlap; self.parallelPartitioningType = parallelPartitioningType
        if self.useC:
            self.meshList.append(EdgeMesh())
            self.meshList[0].generateEdgeMeshFromRectangularGrid(nx,Lx)
            self.cmultilevelMesh = cmeshTools.CMultilevelMesh(self.meshList[0].cmesh,refinementLevels)
            self.buildFromC(self.cmultilevelMesh)
            self.meshList[0].nodeArray[:,0] += x
            self.meshList[0].nodeArray[:,1] += y
            self.meshList[0].nodeArray[:,2] += z
            self.meshList[0].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)
            for l in range(1,refinementLevels):
                self.meshList.append(EdgeMesh())
                self.meshList[l].cmesh = self.cmeshList[l]
                self.meshList[l].buildFromC(self.meshList[l].cmesh)
                self.meshList[l].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)
        else:
            grid=RectangularGrid(nx,ny,nz,Lx,Ly,Lz)
            self.meshList.append(EdgeMesh())
            self.meshList[0].rectangularToEdge(grid)
            self.meshList[0].nodeArray[:,0] += x
            self.meshList[0].nodeArray[:,1] += y
            self.meshList[0].nodeArray[:,2] += z
            self.elementChildren=[]
            print self.meshList[0].meshInfo()
            for l in range(1,refinementLevels):
                self.refine()
                print self.meshList[-1].meshInfo()
    def refine(self):
        self.meshList.append(EdgeMesh())
        childrenDict = self.meshList[-1].refine(self.meshList[-2])
        self.elementChildren.append(childrenDict)
    def computeGeometricInfo(self):
        for m in self.meshList:
            m.computeGeometricInfo()
    def locallyRefine(self,elementTagArray):
        """
        simple local refinement assuming elementTagArray[eN]=1 --> bisect
        """
        if self.useC:
            self.cmeshTools.locallyRefineMultilevelMesh(1,self.cmultilevelMesh,elementTagArray)
            self.buildFromC(self.cmultilevelMesh)
            self.meshList.append(EdgeMesh())
            self.meshList[self.nLevels-1].cmesh = self.cmeshList[self.nLevels-1]
            self.meshList[self.nLevels-1].buildFromC(self.meshList[self.nLevels-1].cmesh)
            self.meshList[self.nLevels-1].partitionMesh(nLayersOfOverlap=self.nLayersOfOverlap,parallelPartitioningType=self.parallelPartitioningType)
        else:
            print """locallyRefine not implemented for self.useC= %s """ % (self.useC)
        #
    #
class MultilevelSimplicialMesh(MultilevelMesh):
    """A wrapper for all the simplicial hierarchical meshes in 1,2, and 3D"""
    def __init__(self,nd,nx,ny=1,nz=1,Lx=1.0,Ly=1.0,Lz=1.0,refinementLevels=1):
        if nd==1:
            MultilevelEdgeMesh.__init__(self,nx,ny,nz,
                                        Lx,Ly,Lz,
                                        refinementLevels)
        elif nd==2:
            MultilevelTriangularMesh.__init__(self,nx,ny,nz,
                                              Lx,Ly,Lz,
                                              refinementLevels)
        elif nd==3:
            MultilevelTetrahedralMesh.__init__(self,nx,ny,nz,
                                               Lz,Ly,Lz,
                                               refineMentLevels)
    def refine(self):
        if nd==1:
            MultilevelEdgeMesh.refine(self)
        elif nd==2:
            MultilevelTrianglularMesh.refine(self)
        elif nd==3:
            MultilevelTetrahedralMesh.refine(self)

## @}

###utility functions for reading meshes from Xdmf
from xml.etree import ElementTree as ET
import tables,os

def findXMLgridElement(xmf,MeshTag='Spatial_Domain',id_in_collection=-1,verbose=0):
    """Try to find the element of the xml tree xmf that holds a uniform
    grid with the name given in MeshTag by searching through Temporal
    Grid Collections and Grid Collections.

    If MeshTag isn't found, uses the first entry in the Domain
    """
    Domain = xmf.getroot()[-1]
    GridCollection = None
    Grid = None
    for collection in Domain:
        if 'Name' in collection.attrib and MeshTag in collection.attrib['Name']:
            GridCollection = collection
            break
    if GridCollection is None:
        GridCollection = Domain[0]
    logEvent("Trying GridCollection.tag= %s" % (GridCollection.tag),4)
    if GridCollection.attrib['GridType'] == 'Collection':
        Grid = GridCollection[-1]
    elif GridCollection.attrib['GridType'] == 'Uniform':
        Grid = GridCollection
    assert Grid.tag == 'Grid'
    assert Grid.attrib['GridType'] == 'Uniform'

    return Grid

def extractPropertiesFromXdmfGridNode(Grid):
    """unpack the Topology, Geometry, NodeMaterials, and ElementMaterials
    nodes from xdmf node for a uniform grid
    """
    #Geometry first
    Topology = None; Geometry  = None; NodeMaterials= None; ElementMaterials = None
    for i,leaf in enumerate(Grid):
        logEvent("Grid leaf %d tag= %s " % (i,leaf.tag),4)
        if leaf.tag == 'Topology':
            Topology = Grid[i]
            logEvent("Topology found in leaf %d " % i,4)
        elif leaf.tag == 'Geometry':
            Geometry = Grid[i]
            logEvent("Geometry found in leaf %d " % i,4)
        elif leaf.tag == 'Attribute' and leaf.attrib['Name'] == 'nodeMaterialTypes':
            NodeMaterials = Grid[i]
            logEvent("NodeMaterials found in leaf %d " % i,4)
        elif leaf.tag == 'Attribute' and leaf.attrib['Name'] == 'elementMaterialTypes':
            ElementMaterials = Grid[i]
            logEvent("ElementMaterials found in leaf %d " % i,4)

    return Topology,Geometry,NodeMaterials,ElementMaterials

def readUniformElementTopologyFromXdmf(elementTopologyName,Topology,hdf5,topologyid2name,topology2nodes):
    """
    Read xmdf element topology information when there are uniform elements in the mesh
    Type of element given by elementTopologyName
    Heavy data stored in hdf5
    topologyid2name -- lookup for number of nodes in a given element type

    returns

    nElements_global  -- the number of elements in the mesh
    nNodes_element    -- number of nodes per element
    elementNodesArray -- element --> node connectivity stored as flattened array accessed using elementNodes_offset
    elementNodes_offset -- offsets into the elementNodesArray storage for element connectivity, 
    element eN nodes are in elementNodesArray[elementNodes_offset[eN]:elementNodes_offset[eN+1]]

    """

    nNodes_element = topology2nodes[elementTopologyName]
    entry = Topology[0].text.split(':')[-1]
    logEvent("Reading  elementNodesArray from %s " % entry,3)

    elementNodesArray = hdf5.get_node(entry).read()
    assert elementNodesArray.shape[1] == nNodes_element
    nElements_global = elementNodesArray.shape[0]
    logEvent("nElements_global,nNodes_element= (%d,%d) " % (nElements_global,nNodes_element),3)

    elementNodes_offset = np.arange(nElements_global*nNodes_element+1,step=nNodes_element,dtype='i')

    return nElements_global, nNodes_element, elementNodesArray, elementNodes_offset

def readMixedElementTopologyFromXdmf(elementTopologyName,Topology,hdf5,topologyid2name,topology2nodes):
    """
    Read xmdf element topology information when there are mixed elements in the mesh
    Heavy data stored in hdf5
    topologyid2name -- lookup for number of nodes in a given element type

    returns

    nElements_global  -- the number of elements in the mesh
    elementNodesArray -- element --> node connectivity stored as flattened 
    array accessed using elementNodes_offset
    elementNodes_offset -- offsets into the elementNodesArray storage for element 
    connectivity, element eN nodes are 
    inelementNodesArray[elementNodes_offset[eN]:elementNodes_offset[eN+1]]
 
    """
    assert elementTopologyName == 'Mixed'

    entry = Topology[0].text.split(':')[-1]
    logEvent("Reading xdmf_topology from %s " % entry,3)

    xdmf_topology = hdf5.get_node(entry).read()
    #build elementNodesArray and offsets now
    nElements_global = 0
    i = 0
    while i < len(xdmf_topology):
        nElements_global += 1
        nNodes_local = topology2nodes[topologyid2name[xdmf_topology[i]]]
        i += nNodes_local+1
    #
    logEvent("Mixed topology found %s elements " % nElements_global,3)
    elementNodes_offset = np.zeros((nElements_global+1,),'i')

    i = 0; eN = 0
    while i < len(xdmf_topology):
        nNodes_local = topology2nodes[topologyid2name[xdmf_topology[i]]]
        elementNodes_offset[eN+1] = elementNodes_offset[eN] + nNodes_local
        eN += 1; i += nNodes_local+1
    elementNodesArray = np.zeros((elementNodes_offset[nElements_global],),'i')
    i = 0; eN = 0
    while i < len(self.xdmf_topology):
        nNodes_local = topology2nodes[topologyid2name[xdmf_topology[i]]]
        elementNodesArray[elementNodes_offset[eN]:elementNodes_offset[eN+1]][:] = xdmf_topology[i+1:i+1+nNodes_local][:]
        eN += 1; i += nNodes_local+1

    return nElements_global, elementNodesArray, elementNodes_offset

def readMeshXdmf(xmf_archive_base,heavy_file_base,MeshTag="Spatial_Domain",hasHDF5=True,verbose=0):
    """Read in a mesh from XDMF, assuming heavy data is in hdf5

    :return: a BasicMeshInfo object with the minimal information read
    
    """    
    # start trying to read an xdmf archive with name xmf_archive_base.xmf
    # assumes heavy_file_base.h5 has heavy data
    # root Element is Xdmf
    #   last child of Xdmf which should be a Domain Element
    #      find child of Domain that is a Temporal Grid Collection with a name containing MeshTag, if None use first collection
    #         last child of Temporal Grid Collection should be a Uniform Grid at final time
    #            Attribute (usually 1) of child is  Topology
    #               set elementTopologyName to Type
    #               if Type != Mixed
    #                 get text attribute and read this entry from  hdf5 file
    #                 set nNodes_element based on Type, nElements_global from leading dimension of elementNodesArray
    #                 create elementNodes_offset from Type and flatten elementNodesArray
    #               else
    #                 get text attribute and read this entry from  hdf5 file to place in into xdmf_topology
    #                 generate elementNodesArray from xdmf_topology, calculating the number of elements using
    #                   walk through xdmf_topology
    #            Attribute (usually 2) of child is Geometry  --> load data into nodeArray
    #                set nNodes_global from nodeArray
    #            If has Attribute nodeMaterials read this from hdf file, else set to default of all zeros
    #            If has Attribute elementMaterialTypes, read this from hdf file, else set to default of all zeros
    assert os.path.isfile(xmf_archive_base+'.xmf')
    assert os.path.isfile(heavy_file_base+'.h5')

    ###information about allowed Xdmf topologies
    #Xdmf cell type id to Name
    topologyid2name = {2:'Polyline',4:'Triangle',5:'Quadrilateral',6:'Tetrahedron',8:'Wedge',9:'Hexahedron',
                       112:'Mixed'} #Mixed isn't actually used 0x070
    #Topology name to number of local nodes
    topology2nodes = {'Polyline':2,'Triangle':3,'Quadrilateral':4,'Tetrahedron':4,'Wedge':6,'Hexahedron':8}

    #for output
    class BasicMeshInfo:
        def __init__(self):
            self.nNodes_global     = None
            self.nodeArray         = None
            self.nodeMaterialTypes = None
            self.nNodes_element    = None
            self.nElements_global  = None
            self.elementTopologyName = None
            self.elementNodesArray = None
            self.elementNodes_offset  = None
            self.elementMaterialTypes = None
            self.nNodes_owned         = None
            self.nElements_owned      = None
        #
    #
    MeshInfo = BasicMeshInfo()

    xmf = ET.parse(xmf_archive_base+'.xmf')
    hdf5= tables.open_file(heavy_file_base+'.h5',mode="r")
    assert hasHDF5

    Grid = findXMLgridElement(xmf,MeshTag,id_in_collection=-1,verbose=verbose)

    Topology,Geometry,NodeMaterials,ElementMaterials = extractPropertiesFromXdmfGridNode(Grid)

    assert Geometry is not None
    entry = Geometry[0].text.split(':')[-1]
    logEvent("Reading nodeArray from %s " % entry,3)

    MeshInfo.nodeArray = hdf5.get_node(entry).read()
    MeshInfo.nNodes_global = MeshInfo.nodeArray.shape[0]

    if NodeMaterials is not None:
        entry = NodeMaterials[0].text.split(':')[-1]
        logEvent("Reading nodeMaterialTypes from %s " % entry,4)
        MeshInfo.nodeMaterialTypes = hdf5.get_node(entry).read()
    else:
        MeshInfo.nodeMaterialTypes = np.zeros((MeshInfo.nNodes_global,),'i')

    assert Topology is not None
    if 'Type' in Topology.attrib:
        MeshInfo.elementTopologyName = Topology.attrib['Type']
    elif 'TopologyType' in Topology.attrib:
        MeshInfo.elementTopologyName = Topology.attrib['TopologyType']
    assert MeshInfo.elementTopologyName is not None

    logEvent("elementTopologyName= %s " % MeshInfo.elementTopologyName,3)
    assert MeshInfo.elementTopologyName in topologyid2name.values()

    if MeshInfo.elementTopologyName != 'Mixed':
        MeshInfo.nElements_global, MeshInfo.nNodes_element, \
            MeshInfo.elementNodesArray, MeshInfo.elementNodes_offset = readUniformElementTopologyFromXdmf(MeshInfo.elementTopologyName,Topology,
                                                                                                          hdf5,topologyid2name,topology2nodes)

    else:
        MeshInfo.nElements_global, MeshInfo.elementNodesArray, \
            MeshInfo.elementNodes_offset = readMixedElementTopologyFromXdmf(MeshInfo.elementTopologyName,Topology,hdf5,topologyid2name,topology2nodes)

    #
    if ElementMaterials is not None:
        entry = ElementMaterials[0].text.split(':')[-1]
        logEvent("Reading elementMaterialTypes from %s " % entry,3)
        MeshInfo.elementMaterialTypes = hdf5.get_node(entry).read()

    else:
        MeshInfo.elementMaterialTypes = np.zeros((MeshInfo.nElements_global,),'i')
    #
    ###only serial for now
    MeshInfo.nNodes_owned = MeshInfo.nNodes_global
    MeshInfo.nElements_owned = MeshInfo.nElements_global
    hdf5.close()

    return MeshInfo
#
def writeHexMesh(mesh_info,hexfile_base,index_base=0):
    """
    Write a hex mesh in Ido's format with base numbering index_base
    HEX
    nNodes_global nElements_global
    x0 y0 z0
    x1 y1 z1
    ...
    xN yN zN
    [n0 n1 n2 n3 n4 n5 n6 n7 mat0]
    [n0 n1 n2 n3 n4 n5 n6 n7 mat1]
    """
    assert mesh_info.elementTopologyName=='Hexahedron'

    header="""HEX
{nNodes_global} {nElements_global}
""".format(nNodes_global=mesh_info.nNodes_global,nElements_global=mesh_info.nElements_global)

    with open(hexfile_base+'.mesh','w') as mout:
        mout.write(header)
        np.savetxt(mout,mesh_info.nodeArray)
        #format the elements, appending element material type
        elems_with_mat = np.append(mesh_info.elementNodesArray,mesh_info.elementMaterialTypes.reshape(mesh_info.nElements_global,1),axis=1)
        elems_with_mat[:,:-1] += index_base
        np.savetxt(mout,elems_with_mat,fmt='%d')


class MultilevelNURBSMesh(MultilevelMesh):
    def __init__(self,
                 nx, ny, nz,
                 x=0.0, y=0.0, z=0.0,
                 px=1, py=1, pz=1,
                 Lx=1.0, Ly=1.0, Lz=1.0,
                 refinementLevels=1,
                 skipInit=False,
                 nLayersOfOverlap=1,
                 parallelPartitioningType=MeshParallelPartitioningTypes.node):
        import cmeshTools
        import Comm
        MultilevelMesh.__init__(self)
        self.useC = True
        self.nLayersOfOverlap = nLayersOfOverlap; self.parallelPartitioningType = parallelPartitioningType
        logEvent("Generating NURBS mesh")
        if not skipInit:
            self.meshList.append(NURBSMesh())
            self.meshList[0].generateNURBSMeshFromRectangularGrid(nx,ny,nz,px,py,pz,Lx,Ly,Lz)
            self.cmultilevelMesh = cmeshTools.CMultilevelMesh(self.meshList[0].cmesh,refinementLevels)
            self.buildFromC(self.cmultilevelMesh)
            self.meshList[0].nodeArray[:,0] += x
            self.meshList[0].nodeArray[:,1] += y
            self.meshList[0].nodeArray[:,2] += z
            self.meshList[0].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)
            for l in range(1,refinementLevels):
                self.meshList.append(NURBSMesh())
                self.meshList[l].cmesh = self.cmeshList[l]
                self.meshList[l].buildFromC(self.cmeshList[l])
                self.meshList[l].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)

    def generateFromExistingCoarseMesh(self,mesh0,refinementLevels,nLayersOfOverlap=1,
                                       parallelPartitioningType=MeshParallelPartitioningTypes.node):
        import cmeshTools
        #blow away or just trust garbage collection
        self.nLayersOfOverlap=nLayersOfOverlap;self.parallelPartitioningType=parallelPartitioningType
        self.meshList = []
        self.elementParents = None
        self.cmultilevelMesh = None

        self.meshList.append(mesh0)
        self.cmultilevelMesh = cmeshTools.CMultilevelMesh(self.meshList[0].cmesh,refinementLevels)
        self.buildFromC(self.cmultilevelMesh)
        self.meshList[0].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)
        for l in range(1,refinementLevels):
            self.meshList.append(NURBSMesh())
            self.meshList[l].cmesh = self.cmeshList[l]
            self.meshList[l].buildFromC(self.meshList[l].cmesh)
            self.meshList[l].partitionMesh(nLayersOfOverlap=nLayersOfOverlap,parallelPartitioningType=parallelPartitioningType)


    def refine(self):
        self.meshList.append(NURBSMesh())
        childrenDict = self.meshList[-1].refine(self.meshList[-2])
        self.elementChildren.append(childrenDict)
    def computeGeometricInfo(self):
        for m in self.meshList:
            m.computeGeometricInfo()





class NURBSMesh(HexahedralMesh):
    """A mesh consisting of NURBS.

    """
    def __init__(self):
        HexahedralMesh.__init__(self)

    def generateHexahedralMeshFromRectangularGrid(self,nx,ny,nz,Lx,Ly,Lz):
        generateNURBSMeshFromRectangularGrid(self,nx,ny,nz,1,1,1,Lx,Ly,Lz)

    def generateNURBSMeshFromRectangularGrid(self,nx,ny,nz,px,py,pz,Lx,Ly,Lz):
        import cmeshTools
        self.cmesh = cmeshTools.CMesh()
        cmeshTools.generateNURBSMeshFromRectangularGrid(nx,ny,nz,px,py,pz,Lx,Ly,Lz,self.cmesh)
        cmeshTools.allocateGeometricInfo_NURBS(self.cmesh)
        cmeshTools.computeGeometricInfo_NURBS(self.cmesh)
        self.buildFromC(self.cmesh)


def distance(a, b):
    norm = np.linalg.norm
    return norm(b - a)

def triangleVerticesToNormals(elementVertices):
    """
    Given a set of vertices to a triangle, return normals and a point corresponding to each normal
    """
    norm = np.linalg.norm

    elementVertices = np.asarray(elementVertices)

    if norm(elementVertices[:,2]) > 0:
        raise ValueError("Expected triangles in 2D plane, got something else")
    sets = ((0, 1), (0, 2), (1, 2))
    outs = (2, 1, 0)
    faces = []
    rotate = np.asarray(((0., -1., 0.),
                         (1., 0., 0.),
                         (0., 0., 0.)))

    for set, out in zip(sets, outs):
        vertices = elementVertices[[set]]
        ab = vertices[1] - vertices[0]
        v_out = vertices[0] - elementVertices[out]
        normal = rotate.dot(ab)
        # normal should point *away* from remaining point
        if normal.dot(v_out) < 0:
            normal = -1*normal
        faces.append((normal, vertices[0]))
    return faces

def tetrahedronVerticesToNormals(elementVertices):
    """
    Given a set of vertices to a tetrahedron, return normals and a point corresponding to each normal
    """

    elementVertices = np.asarray(elementVertices)

    sets = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
    outs = (3, 2, 1, 0)

    faces = []

    for set, out in zip(sets, outs):
        vertices = elementVertices[[set]]
        ab = vertices[1] - vertices[0]
        ac = vertices[2] - vertices[0]
        normal = np.cross(ab, ac)
        v_out = vertices[0] - elementVertices[out]
        # normal should point *away* from remaining point
        if normal.dot(v_out) < 0:
            normal = -1*normal
        faces.append((normal, vertices[0]))
    return faces


def intersectPoints(line, points):
    """
    Given a line segment (defined as two points), identify all points that the line segment intersects.

    This hasn't been vectorized.
    """

    a, b = line
    a = np.asarray(a)
    b = np.asarray(b)
    distanceAB = distance(a, b)

    def onAB(p):
        p = np.asarray(p)
        eps = 2*np.max((np.max(np.spacing(a)), np.max(np.spacing(b)), np.max(np.spacing(p))))
        distancePA = distance(a, p)
        distancePB = distance(p, b)
        return p if abs(distancePA + distancePB - distanceAB) < eps else None

    return [onAB(p) for p in points]


def intersectEdges(line, edges):
    """
    Given a line segment (defined as two points), identify the locations of its intersections with all
    given edges (defined as line segments).  If the line and an edge overlap, the *furthest* point
    along the line (closest to the second point) that is still on each edge is returned.

    This hasn't been vectorized.
    """
    norm = np.linalg.norm

    def intersectEdge(line, edge):

        line = np.asarray(line)
        edge = np.asarray(edge)
        a, b = line
        c, d = edge
        v_l = b - a
        v_e = d - c

        vl_cross_ve = np.cross(v_l, v_e)
        mag_vl_cross_ve = norm(vl_cross_ve)

        if mag_vl_cross_ve == 0:
            # lines are parallel, check for overlap
            intersects = intersectPoints(line, edge) + intersectPoints(edge, line)
            # test for an intersect in intersectPoints
            intersect = next((i for i in intersects if i is not None), None)
            if intersect is not None:
                # farthest endpoint is a, so start from there
                closest_endpoint = a
                closest_distance = distance(closest_endpoint, b)
                # could reuse iterator from above, but it's confusing enough as it is :)
                for intersect in intersects:
                    if intersect is None:
                        continue
                    intersect_distance = distance(intersect, b)
                    if intersect_distance < closest_distance:
                        closest_endpoint = intersect
                        closest_distance = intersect_distance
                return closest_endpoint
            else:
                return None

        # lines are not parallel, check for intersection
        vl_cross_ve = np.cross(v_l, v_e)

        # if v_l and v_e intersect, then there is an x that satisfies
        x_vl_cross_ve = np.cross((c - a), v_e)

        # but the two above vectors must be parallel
        if norm(np.cross(vl_cross_ve, x_vl_cross_ve)) > 1e-8:
            return None

        # two lines are parallel, solve for x
        x = norm(x_vl_cross_ve)/norm(vl_cross_ve)

        intersect = a + x*(b-a)

        # and verify intersection is on the line
        points = intersectPoints(line, [intersect])
        assert(len(points) == 1)
        return points[0]

    return [intersectEdge(line, edge) for edge in edges]


def intersectPolyhedron(line, polyhedron):
    """
    Given a line (defined as two points), identify the locations that it enters and exits the
    polyhedron (defined as a collection of half-planes in three-space in normal, vertex form)

    If the facets of the polyhedron are in edge form, the normal can be computed by taking the cross product of any
    two non-parallel edges of the facet (in three-space).  Any vertex of the facet will work.

    Implementation of algorithm described here: http://geomalgorithms.com/a13-_intersect-4.html

    This hasn't been vectorized.
    """

    a, b = line
    a, b = np.asarray(a), np.asarray(b)

    if distance(a, b) == 0:
        raise ValueError("Line segment must not have length 0")

    v_l = b - a
    t_e = 0  # location along line entering polyhedron (initial value 0)
    t_l = 1  # location along line leaving polyhedron (initial value 1)

    for plane in polyhedron:
        n, v = plane
        n, v = np.asarray(n), np.asarray(v)
        ndotba = -n.dot(a - v)
        d = n.dot(v_l)
        if d == 0:
            # the line segment is parallel to this face
            if ndotba < 0:
                # the line is outside the face
                return None
            else:
                # the line is in or on the face, ignore this face
                continue
        t = ndotba / float(d)
        if d < 0:
            # segment is entering polyhedron across this facet
            t_e = max(t_e, t)
            if t_e > t_l:
                # segment enters polyhedron after leaving, no intersection
                return None
        else:
            # segment is exiting polyhedron across this facet
            t_l = min(t_l, t)
            if t_l < t_e:
                # segment exits polyhedron before entering, no intersection
                return None

    assert(t_e <= t_l)

    return [a + t_e*v_l, a + t_l*v_l]


def getMeshIntersections(mesh, toPolyhedron, endpoints):
    """
    Return all intersections between a line segment and a Proteus mesh


    :param mesh - a Proteus mesh
    :param toPolyhedron - a method for converting Proteus element vertices to polyhedra in normal/point form
    :param endpoints - a pair of points in 3-space defining the line segment

    :return a list of pairs of intersections through the mesh
    """

    intersections = set()
    for element in mesh.elementNodesArray:
        # map nodes to physical vertices
        elementVertices = mesh.nodeArray[element]
        # get plane normals
        polyhedron = toPolyhedron(elementVertices)
        elementIntersections = intersectPolyhedron(endpoints, polyhedron)
        if elementIntersections:
            if np.array_equal(elementIntersections[0], elementIntersections[1]):
                continue
            intersections.update(((tuple(elementIntersections[0]), tuple(elementIntersections[1])),),)
    return intersections

def runTetgen(polyfile,
              baseFlags="Yp",
              name = ""):
    """
    Generate tetgen files from a polyfile.

    Arguments
    ---------
    polyfile : str
        Filename with appropriate data for tengen.
    baseFlags : str
        Standard Tetgen options for generation
    name : str
        

    """
    from subprocess import check_call
    tetcmd = "tetgen - %s %s.poly" % (baseFlags, polyfile)
    
    check_call(tetcmd,shell=True)
    
    logEvent("Done running tetgen")
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

def genMeshWithTetgen(polyfile,
                      nbase=1):
   """
   Generate a mesh from a set of tetgen files.

   Arguments
   ---------
   polyfile : str
       Filename base for tetgen files
   nbase : int

   Returns
   --------
   mesh : :class:`proteus.MeshTools.TetrahedralMesh`
       Simplex mesh
   """
   elefile = "%s.ele" % polyfile
   nodefile = "%s.node" % polyfile
   facefile = "%s.face" % polyfile
   edgefile = "%s.edge" % polyfile
   assert os.path.exists(elefile), "no .ele file"
   assert os.path.exists(nodefile), "no  .node file"
   assert os.path.exists(facefile), "no .face file"
   mesh=TetrahedralMesh()
   mesh.generateFromTetgenFiles(polyfile,
                                base=nbase)
   return mesh

class MeshOptions:
    """
    Mesh options for the domain

    Parameters
    ----------
    domain: proteus.Domain
    """
    def __init__(self, domain):
        self.Domain = domain
        self.he = 1.
        self.use_gmsh = False
        self.genMesh = True
        self.outputFiles_name = 'mesh'
        self.outputFiles = {'poly': True,     
                            'ply': False,        
                            'asymptote': False,
                            'geo': False}
        self.restrictFineSolutionToAllMeshes = False
        self.parallelPartitioningType = MeshParallelPartitioningTypes.node
        self.nLayersOfOverlapForParallel = 1
        self.triangleOptions = "q30DenA" # defined when setTriangleOptions called
        self.nLevels = 1
        if domain is not None:
            self.nd = domain.nd
            if self.nd == 2:
                self.triangle_string = 'VApq30Dena'
            if self.nd == 3:
                self.triangle_string = 'VApq1.35q12feena'
        else:
            self.triangle_string = None

    def setElementSize(self, he):
        """
        Sets element size for uniform mesh.

        Parameters
        ----------
        he: float
            mesh characteristic element size
        """
        self.he = he

    def setParallelPartitioningType(self, partitioning_type='node', layers_overlap=0):
        """
        Changes parallel partitioning type

        Parameters
        ----------
        partitioning_type: Optional[str, int]
            parallel partitioning type (default: 'node' (1))
        layers: int
            layers of overlap for paralllel (default: 0)
        """
        if partitioning_type == 'element' or partitioning_type == 0:
            self.parallelPartitioningType = MeshParallelPartitioningTypes.element
        if partitioning_type == 'node' or partitioning_type == 1:
            self.parallelPartitioningType = MeshParallelPartitioningTypes.node
        self.nLayersOfOverlapForParallel = layers_overlap

    def setTriangleOptions(self, triangle_options=None):
        """
        Sets the trangle options

        Parameters
        ----------
        triangle_options: Optional[str]
            string for triangle options. If not passed, it will be
            set with triangle_string attribute and 'he' value, with
            default for 2D: he**2/2; default for 3D: he**3/6
        """
        if triangle_options is not None:
            self.triangleOptions = triangle_options
        else:
            assert self.he is not None, 'Element size (he) must be set before setting triangle options'
            assert self.triangle_string is not None, 'triangle_string must be set before setting triangle options'
            if self.nd == 2:
                self.triangleOptions = self.triangle_string + '%8.8f' \
                                       % (self.he**2/2.,)
            elif self.nd == 3:
                self.triangleOptions = self.triangle_string + '%21.16e' \
                                       % (self.he**3/6.,)

    def setMeshGenerator(self, generator):
        """
        Indicates mesh generator to use

        Parameters
        ----------
        generator: str
            options: 'gmsh', 'triangle', 'tetgen'

        (!) Only has an effect when setting to 'gmsh' in current 
        implementation (triangle is default for 2D, tetgen for 3D)
        """
        generators = ['gmsh', 'triangle', 'tetgen']
        assert generator in generators, 'Unknown mesh generator'
        if generator == 'gmsh':
            self.use_gmsh = True
        else:
            self.use_gmsh = False

    def setOutputFiles(self, name='mesh', poly=True, ply=False, asymptote=False, geo=False):
        """
        Output files to be created 

        Parameters
        ----------
        name: Optional[str]
            name of the mesh files (prefix) (default: 'mesh')
        poly: Optional[bool]
            create a poly file
        ply: Optional[bool]
            create a ply file
        asymptote: Optional[bool]
            create an asymptote file
        geo:
            create a geofile
        """
        self.outputFiles_name = name
        self.outputFiles['poly'] = poly
        self.outputFiles['ply'] = ply
        self.outputFiles['asymptote'] = asymptote
        self.outputFiles['geo'] = geo


def msh2simplex(fileprefix, nd):
    """
    Converts a .msh file (Gmsh) to .ele .edge .node files (triangle).
    (!) Works only with triangle elements in 2D and tetrahedral elements in 3D.

    Parameters
    ----------
    fileprefix: str
        prefix of the .msh file (e.g. 'mesh' if file called 'mesh.msh')

    """
    assert nd == 2 or nd == 3, 'nd must be 2 or 3'
    mshfile = open(fileprefix+'.msh', 'r')
    nodes = []
    edges_msh = []
    triangles = []
    tetrahedra = []
    tetrahedron_nb = 0
    triangle_nb = 0
    edge_nb = 0
    switch = None
    switch_count = -1
    logEvent('msh2simplex: getting nodes and elements')
    for i, line in enumerate(mshfile):
        if 'Nodes' in line:
            switch = 'nodes'
            switch_count = -1
        if 'Elements' in line:
            switch = 'elements'
            switch_count = -1
        if switch == 'nodes' and switch_count >= 0:
            words = line.split()
            if switch_count == 0:
                node_nb = int(words[0])
            else:
                nid = int(words[0])
                if nd == 2:
                    x, y, z = float(words[1]), float(words[2]), 0
                elif nd == 3:
                    x, y, z = float(words[1]), float(words[2]), float(words[3])
                nodes += [[nid, x, y, z, 0]]
        if switch == 'elements' and switch_count >= 0:
            words = line.split()
            if switch_count == 0:
                el_nb = int(words[0])
            else:
                el_id = int(words[0])
                el_type = int(words[1])
                nb_tags = int(words[2])
                if nb_tags == 2:
                    flag = int(words[3])
                else:
                    flag = 0
                s = 3+nb_tags # starting index on words for element info
                if el_type == 1: # segment
                    edge_nb += 1
                    edges_msh += [[edge_nb, int(words[s]), int(words[s+1]), flag]]
                elif el_type == 2: # triangle
                    triangle_nb += 1
                    triangles += [[triangle_nb, int(words[s]), int(words[s+1]), int(words[s+2]), flag]]
                    # update nodes flags
                    if nd == 3:
                        for i in range(3):
                            if nodes[int(words[s+i])-1][4] == 0:
                                nodes[int(words[s+i])-1][4] = flag
                elif el_type == 4: # tetrahedron 
                    tetrahedron_nb += 1
                    tetrahedra += [[tetrahedron_nb, int(words[s]), int(words[s+1]), int(words[s+2]), int(words[s+3]), flag]]
                elif el_type == 15: # node
                    nodes[el_id-1][4] = flag
        switch_count += 1
    mshfile.close()

    # construct ALL edges with flags and add flags to nodes
    edges_dict = {}
    triangles = np.array(triangles)
    edge_nb = 0
    edges = []

    logEvent('msh2simplex: constructing edges')
    for triangle in triangles[:,1:4]:  # take only vertices index
        for i in range(len(triangle)):
            edge = Edge(edgeNumber=edge_nb, nodes=[triangle[i-1], triangle[i]])
            edge_exist = bool(edges_dict.get(edge.nodes))
            if not edge_exist:
                edge_nb += 1
                edges_dict[edge.nodes] = edge
                edges += [[edge_nb, edge.nodes[0], edge.nodes[1], 0]]
    logEvent('msh2simplex: updating edges and nodes flags')
    edges = np.array(edges)
    for edge in edges_msh:
        edge_nodes = [edge[1], edge[2]]
        edge_nodes.sort()
        edge_nodes = tuple(edge_nodes)
        edge_class = edges_dict.get(edge_nodes)
        edges[edge_class.N, 3] = edge[3]
        # ! edge nodes are indexed from 1 with gmsh
        if nodes[edge[1]-1][-1] == 0:  # update node flags
            nodes[edge[1]-1][-1] = edge[3]
        if nodes[edge[2]-1][-1] == 0:  # update node flags
            nodes[edge[1]-1][-1] = edge[3]

    if nd == 2:
        logEvent('msh2simplex: writing .node .ele .edge files')
    elif nd == 3:
        logEvent('msh2simplex: writing .node .ele .edge .face files')
    header = '{0:d} {1:d} 0 1'.format(node_nb, nd)

    if nd == 2:
        nodes = np.array(nodes)
        nodes = np.delete(nodes, 3, 1)
        fmt = ['%d', '%f', '%f', '%d']
    elif nd == 3:
       fmt = ['%d', '%f', '%f', '%f', '%d']
    np.savetxt(fileprefix+'.node', nodes, fmt=fmt, header=header, comments='')

    header = '{0:d} 1'.format(edge_nb)
    np.savetxt(fileprefix+'.edge', edges, fmt='%d', header=header, comments='')

    if nd == 2:
        header = '{0:d} 3 1'.format(triangle_nb)
        np.savetxt(fileprefix+'.ele', triangles, fmt='%d', header=header, comments='')
    elif nd == 3:
        header = '{0:d} 3 1'.format(triangle_nb)
        np.savetxt(fileprefix+'.face', triangles, fmt='%d', header=header, comments='')
        header = '{0:d} 4 1'.format(tetrahedron_nb)
        np.savetxt(fileprefix+'.ele', tetrahedra, fmt='%d', header=header, comments='')

    logEvent('msh2simplex: finished converting .msh to simplex files')
