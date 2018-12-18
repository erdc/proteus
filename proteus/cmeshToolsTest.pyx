# A type of -*- python -*- file
import cython
cimport cython
import numpy as np
cimport numpy as np
from libcpp cimport bool
from mpi4py import MPI
cimport mesh as cppm
from proteus import Comm

cdef class CMeshLink:
    def __init__(self):
        self.mesh = cppm.Mesh()

# @cython.auto_pickle(True)
cdef class CMesh:
    def __init__(self):
        self.meshlink = CMeshLink()
        cppm.initializeMesh(self.meshlink.mesh)

    def buildPythonMeshInterface(self):
        cdef int dim1
        self.nElements_global = self.meshlink.mesh.nElements_global
        self.nNodes_global = self.meshlink.mesh.nNodes_global
        self.nNodes_element = self.meshlink.mesh.nNodes_element
        self.nNodes_elementBoundary = self.meshlink.mesh.nNodes_elementBoundary
        self.nElementBoundaries_element = self.meshlink.mesh.nElementBoundaries_element
        self.nElementBoundaries_global = self.meshlink.mesh.nElementBoundaries_global
        self.nInteriorElementBoundaries_global = self.meshlink.mesh.nInteriorElementBoundaries_global
        self.nExteriorElementBoundaries_global = self.meshlink.mesh.nExteriorElementBoundaries_global
        self.max_nElements_node = self.meshlink.mesh.max_nElements_node
        self.nEdges_global = self.meshlink.mesh.nEdges_global
        self.max_nNodeNeighbors_node = self.meshlink.mesh.max_nNodeNeighbors_node
        self.elementNodesArray = np.asarray(<int[:self.meshlink.mesh.nElements_global, :self.meshlink.mesh.nNodes_element]> self.meshlink.mesh.elementNodesArray)
        if self.meshlink.mesh.nodeElementOffsets != NULL:
            self.nodeElementsArray = np.asarray(<int[:self.meshlink.mesh.nNodes_global]> self.meshlink.mesh.nodeElementsArray)
        else:
            self.nodeElementsArray = np.empty(0, dtype=np.int32)
        self.nodeElementOffsets = np.asarray(<int[:self.meshlink.mesh.nNodes_global+1]> self.meshlink.mesh.nodeElementOffsets)
        self.elementNeighborsArray = np.asarray(<int[:self.meshlink.mesh.nElements_global, :self.meshlink.mesh.nElementBoundaries_element]> self.meshlink.mesh.elementNeighborsArray)
        self.elementBoundariesArray = np.asarray(<int[:self.meshlink.mesh.nElements_global, :self.meshlink.mesh.nElementBoundaries_element]> self.meshlink.mesh.elementBoundariesArray)
        self.elementBoundaryNodesArray = np.asarray(<int[:self.meshlink.mesh.nElementBoundaries_global, :self.meshlink.mesh.nNodes_elementBoundary]> self.meshlink.mesh.elementBoundaryNodesArray)
        self.elementBoundaryElementsArray = np.asarray(<int[:self.meshlink.mesh.nElementBoundaries_global, :2]> self.meshlink.mesh.elementBoundaryElementsArray)
        self.elementBoundaryLocalElementBoundariesArray = np.asarray(<int[:self.meshlink.mesh.nElementBoundaries_global, :2]> self.meshlink.mesh.elementBoundaryLocalElementBoundariesArray)
        self.interiorElementBoundariesArray = np.asarray(<int[:self.meshlink.mesh.nInteriorElementBoundaries_global]> self.meshlink.mesh.interiorElementBoundariesArray)
        self.exteriorElementBoundariesArray = np.asarray(<int[:self.meshlink.mesh.nExteriorElementBoundaries_global]> self.meshlink.mesh.exteriorElementBoundariesArray)
        self.edgeNodesArray = np.asarray(<int[:self.meshlink.mesh.nEdges_global, :2]> self.meshlink.mesh.edgeNodesArray)
        if self.meshlink.mesh.nodeStarOffsets != NULL:
            self.nodeStarArray = np.asarray(<int[:self.meshlink.mesh.nodeStarOffsets[self.meshlink.mesh.nNodes_global]]> self.meshlink.mesh.nodeStarArray)
        else:
            self.nodeStarArray = np.empty(0, dtype=np.int32)
        self.nodeStarOffsets = np.asarray(<int[:self.meshlink.mesh.nNodes_global+1]> self.meshlink.mesh.nodeStarOffsets)
        self.elementMaterialTypes = np.asarray(<int[:self.meshlink.mesh.nElements_global]> self.meshlink.mesh.elementMaterialTypes)
        self.elementBoundaryMaterialTypes = np.asarray(<int[:self.meshlink.mesh.nElementBoundaries_global]> self.meshlink.mesh.elementBoundaryMaterialTypes)
        self.nodeMaterialTypes = np.asarray(<int[:self.meshlink.mesh.nNodes_global]> self.meshlink.mesh.nodeMaterialTypes)
        self.nodeArray = np.asarray(<double[:self.meshlink.mesh.nNodes_global, :3]> self.meshlink.mesh.nodeArray)
        self.nx = self.meshlink.mesh.nx
        self.ny = self.meshlink.mesh.ny
        self.nz = self.meshlink.mesh.nz
        self.px = self.meshlink.mesh.px
        self.py = self.meshlink.mesh.py
        self.pz = self.meshlink.mesh.pz
        if self.meshlink.mesh.elementIJK != NULL:
            self.elementIJK = np.asarray(<int[:self.meshlink.mesh.nElements_global*3]> self.meshlink.mesh.elementIJK)
        else:
            self.elementIJK = np.empty(0, dtype=np.int32)
        if self.meshlink.mesh.weights != NULL:
            self.weights = np.asarray(<double[:self.meshlink.mesh.nElements_global]> self.meshlink.mesh.weights)
        else:
            self.weights = np.empty(0)
        if self.meshlink.mesh.U_KNOT != NULL:
            self.U_KNOT = np.asarray(<double[:self.meshlink.mesh.nx+self.meshlink.mesh.px+1]> self.meshlink.mesh.U_KNOT)
        else:
            self.U_KNOT = np.empty(0)
        if self.meshlink.mesh.V_KNOT != NULL:
            self.V_KNOT = np.asarray(<double[:self.meshlink.mesh.ny+self.meshlink.mesh.py+1]> self.meshlink.mesh.V_KNOT)
        else:
            self.V_KNOT = np.empty(0)
        if self.meshlink.mesh.W_KNOT != NULL:
            self.W_KNOT = np.asarray(<double[:self.meshlink.mesh.nz+self.meshlink.mesh.pz+1]> self.meshlink.mesh.W_KNOT)
        else:
            self.W_KNOT = np.empty(0)
        self.elementDiametersArray = np.asarray(<double[:self.meshlink.mesh.nElements_global]> self.meshlink.mesh.elementDiametersArray)
        self.elementInnerDiametersArray = np.asarray(<double[:self.meshlink.mesh.nElements_global]> self.meshlink.mesh.elementInnerDiametersArray)
        self.elementBoundaryDiametersArray = np.asarray(<double[:self.meshlink.mesh.nElementBoundaries_global]> self.meshlink.mesh.elementBoundaryDiametersArray)
        self.elementBarycentersArray = np.asarray(<double[:self.meshlink.mesh.nElements_global, :3]> self.meshlink.mesh.elementBarycentersArray)
        self.elementBoundaryBarycentersArray = np.asarray(<double[:self.meshlink.mesh.nElementBoundaries_global, :3]> self.meshlink.mesh.elementBoundaryBarycentersArray)
        self.nodeDiametersArray = np.asarray(<double[:self.meshlink.mesh.nNodes_global]> self.meshlink.mesh.nodeDiametersArray)
        self.nodeSupportArray = np.asarray(<double[:self.meshlink.mesh.nNodes_global]> self.meshlink.mesh.nodeSupportArray)
        self.h = self.meshlink.mesh.h
        self.hMin = self.meshlink.mesh.hMin
        self.sigmaMax = self.meshlink.mesh.sigmaMax
        self.volume = self.meshlink.mesh.volume
    
def buildPythonMeshInterface(cmesh):
    """
    function to be conform to old module, and to calls from MeshTools
    """
    cmesh.buildPythonMeshInterface()
    return (cmesh.nElements_global,
            cmesh.nNodes_global,
            cmesh.nNodes_element,
            cmesh.nNodes_elementBoundary,
            cmesh.nElementBoundaries_element,
            cmesh.nElementBoundaries_global,
            cmesh.nInteriorElementBoundaries_global,
            cmesh.nExteriorElementBoundaries_global,
            cmesh.max_nElements_node,
            cmesh.nEdges_global,
            cmesh.max_nNodeNeighbors_node,
            cmesh.elementNodesArray,
            cmesh.nodeElementsArray,
            cmesh.nodeElementOffsets,
            cmesh.elementNeighborsArray,
            cmesh.elementBoundariesArray,
            cmesh.elementBoundaryNodesArray,
            cmesh.elementBoundaryElementsArray,
            cmesh.elementBoundaryLocalElementBoundariesArray,
            cmesh.interiorElementBoundariesArray,
            cmesh.exteriorElementBoundariesArray,
            cmesh.edgeNodesArray,
            cmesh.nodeStarArray,
            cmesh.nodeStarOffsets,
            cmesh.elementMaterialTypes,
            cmesh.elementBoundaryMaterialTypes,
            cmesh.nodeMaterialTypes,
            cmesh.nodeArray,
            cmesh.nx,
            cmesh.ny,
            cmesh.nz,
            cmesh.px,
            cmesh.py,
            cmesh.pz,
            cmesh.elementIJK,
            cmesh.weights,
            cmesh.U_KNOT,
            cmesh.V_KNOT,
            cmesh.W_KNOT,
            cmesh.elementDiametersArray,
            cmesh.elementInnerDiametersArray,
            cmesh.elementBoundaryDiametersArray,
            cmesh.elementBarycentersArray,
            cmesh.elementBoundaryBarycentersArray,
            cmesh.nodeDiametersArray,
            cmesh.nodeSupportArray,
            cmesh.h,
            cmesh.hMin,
            cmesh.sigmaMax,
            cmesh.volume)

cdef CMesh CMesh_FromMesh(cppm.Mesh mesh):
    CMeshnew = CMesh()
    CMeshnew.meshlink.mesh = mesh
    return CMeshnew

cdef class CMultilevelMeshLink:
    cdef cppm.MultilevelMesh multilevelMesh
    def __init__(self):
        self.multilevelMesh = cppm.MultilevelMesh()
   

cdef class CMultilevelMesh:
    cdef CMultilevelMeshLink multilevelMeshLink
    cdef public:
        int nLevels
        list cmeshList
        list elementParentsArrayList
        list elementChildrenArrayList
        list elementChildrenOffsetsList
    def __init__(self, CMesh cmesh, nLevels):
        self.multilevelMeshLink = CMultilevelMeshLink()
        if cmesh.meshlink.mesh.nNodes_element == 2:
            for i in range(1, nLevels):
                cppm.globallyRefineEdgeMesh(nLevels,
                                            cmesh.meshlink.mesh,
                                            self.multilevelMeshLink.multilevelMesh,
                                            False);
                cppm.allocateGeometricInfo_edge(self.multilevelMeshLink.multilevelMesh.meshArray[i]);
                cppm.computeGeometricInfo_edge(self.multilevelMeshLink.multilevelMesh.meshArray[i]);
                cppm.assignElementBoundaryMaterialTypesFromParent(self.multilevelMeshLink.multilevelMesh.meshArray[i-1],
                                                             self.multilevelMeshLink.multilevelMesh.meshArray[i],
                                                             self.multilevelMeshLink.multilevelMesh.elementParentsArray[i],
                                                             1);
        elif cmesh.meshlink.mesh.nNodes_element == 3:
            cppm.globallyRefineTriangularMesh(nLevels,
                                              cmesh.meshlink.mesh,
                                              self.multilevelMeshLink.multilevelMesh,
                                              False);
            for i in range(1, nLevels):
                cppm.constructElementBoundaryElementsArray_triangle(self.multilevelMeshLink.multilevelMesh.meshArray[i]);
                cppm.allocateGeometricInfo_triangle(self.multilevelMeshLink.multilevelMesh.meshArray[i]);
                cppm.computeGeometricInfo_triangle(self.multilevelMeshLink.multilevelMesh.meshArray[i]);
                cppm.assignElementBoundaryMaterialTypesFromParent(self.multilevelMeshLink.multilevelMesh.meshArray[i-1],
						                                     self.multilevelMeshLink.multilevelMesh.meshArray[i],
						                                     self.multilevelMeshLink.multilevelMesh.elementParentsArray[i],
						                                     2);
        elif cmesh.meshlink.mesh.nNodes_element == 4 and cmesh.meshlink.mesh.nNodes_elementBoundary == 2:
            cppm.globallyRefineQuadrilateralMesh(nLevels,
                                                 cmesh.meshlink.mesh,
                                                 self.multilevelMeshLink.multilevelMesh,
                                                 False);
            for i in range(1, nLevels):
                cppm.constructElementBoundaryElementsArray_quadrilateral(self.multilevelMeshLink.multilevelMesh.meshArray[i]);
                cppm.allocateGeometricInfo_quadrilateral(self.multilevelMeshLink.multilevelMesh.meshArray[i]);
                cppm.computeGeometricInfo_quadrilateral(self.multilevelMeshLink.multilevelMesh.meshArray[i]);
                cppm.assignElementBoundaryMaterialTypesFromParent(self.multilevelMeshLink.multilevelMesh.meshArray[i-1],
						                                          self.multilevelMeshLink.multilevelMesh.meshArray[i],
						                                          self.multilevelMeshLink.multilevelMesh.elementParentsArray[i],
						                                          2);
        elif cmesh.meshlink.mesh.nNodes_element == 4 and cmesh.meshlink.mesh.nNodes_elementBoundary == 3:
            cppm.globallyRefineTetrahedralMesh(nLevels,
                                               cmesh.meshlink.mesh,
                                               self.multilevelMeshLink.multilevelMesh,
                                               False);
            for i in range(1, nLevels):
                cppm.constructElementBoundaryElementsArray_tetrahedron(self.multilevelMeshLink.multilevelMesh.meshArray[i]);
                cppm.allocateGeometricInfo_tetrahedron(self.multilevelMeshLink.multilevelMesh.meshArray[i]);
                cppm.computeGeometricInfo_tetrahedron(self.multilevelMeshLink.multilevelMesh.meshArray[i]);
                cppm.assignElementBoundaryMaterialTypesFromParent(self.multilevelMeshLink.multilevelMesh.meshArray[i-1],
						                                     self.multilevelMeshLink.multilevelMesh.meshArray[i],
						                                     self.multilevelMeshLink.multilevelMesh.elementParentsArray[i],
						                                     3);
        elif cmesh.meshlink.mesh.nNodes_element == 8 and cmesh.meshlink.mesh.nNodes_elementBoundary == 4:
            cppm.globallyRefineHexahedralMesh(nLevels,
                                         cmesh.meshlink.mesh,
                                         self.multilevelMeshLink.multilevelMesh,
                                         False);
            for i in range(1, nLevels):
                cppm.constructElementBoundaryElementsArray_hexahedron(self.multilevelMeshLink.multilevelMesh.
                                                                      meshArray[i]);
                cppm.allocateGeometricInfo_hexahedron(self.multilevelMeshLink.multilevelMesh.meshArray[i]);
                cppm.computeGeometricInfo_hexahedron(self.multilevelMeshLink.multilevelMesh.meshArray[i]);
                cppm.assignElementBoundaryMaterialTypesFromParent(self.multilevelMeshLink.multilevelMesh.meshArray[i-1],
						                                          self.multilevelMeshLink.multilevelMesh.meshArray[i],
						                                          self.multilevelMeshLink.multilevelMesh.elementParentsArray[i],
						                                          4);
        else:
            assert nLevels == 1, 'wrong nLevels'
            cppm.globallyRefineTetrahedralMesh(nLevels,
                                               cmesh.meshlink.mesh,
                                               self.multilevelMeshLink.multilevelMesh,
                                               False);

    def buildPythonMultilevelMeshInterface(self):
        cdef int dim
        self.cmeshList = [CMesh_FromMesh(self.multilevelMeshLink.multilevelMesh.meshArray[0])]
        self.elementParentsArrayList = [np.zeros(0)]
        self.elementChildrenArrayList = []
        self.elementChildrenOffsetsList = []
        for n in range(1, self.multilevelMeshLink.multilevelMesh.nLevels):
            self.cmeshList += [CMesh_FromMesh(self.multilevelMeshLink.multilevelMesh.meshArray[n])]
            dim = self.multilevelMeshLink.multilevelMesh.meshArray[n].nElements_global
            self.elementParentsList += [np.asarray(<int[:dim]> self.multilevelMeshLink.multilevelMesh.elementParentsArray[n])]
            dim = self.multilevelMeshLink.multilevelMesh.elementChildrenOffsets[n-1][self.multilevelMeshLink.multilevelMesh.meshArray[n-1].nElements_global]
            self.elementChildrenList += [np.asarray(<int[:dim]> self.multilevelMeshLink.multilevelMesh.elementChildrenArray[n-1])]
            dim = self.cmultilevelMesh.meshArray[n-1].nElements_global+1;
            self.elementChildrenOffsetsList += [np.asarray(<int[:dim]> self.multilevelMeshLink.multilevelMesh.elementChildrenOffsets[n-1])]
        print(self.cmeshList)

def buildPythonMultilevelMeshInterface(cmultilevelMesh):
    cmultilevelMesh.buildPythonMultilevelMeshInterface()
    return (cmultilevelMesh.nLevels,
            cmultilevelMesh.cmeshList,
            cmultilevelMesh.elementParentsArrayList,
            cmultilevelMesh.elementChildrenArrayList,
            cmultilevelMesh.elementChildrenOffsetsList)

cpdef void generateTetrahedralMeshFromRectangularGrid(int nx,
                                                      int ny,
                                                      int nz,
                                                      double Lx,
                                                      double Ly,
                                                      double Lz,
                                                      CMesh cmesh):
    cppm.regularHexahedralToTetrahedralMeshElements(nx,ny,nz,cmesh.meshlink.mesh);
    cppm.regularHexahedralToTetrahedralMeshNodes(nx,ny,nz,Lx,Ly,Lz,cmesh.meshlink.mesh);
    cppm.constructElementBoundaryElementsArray_tetrahedron(cmesh.meshlink.mesh);
    cppm.regularHexahedralToTetrahedralElementBoundaryMaterials(Lx,Ly,Lz,cmesh.meshlink.mesh);

cpdef void cmeshToolsComputeGeometricInfo_tetrahedron(CMesh cmesh):
    cppm.computeGeometricInfo_tetrahedron(cmesh.meshlink.mesh)

# cdef static void cmeshToolsLocallyRefineMultilevelMesh(CMultilevelMesh cmesh,



cpdef void generateFromTriangleFiles(CMesh cmesh,
                                    const char* filebase,
                                    int base):

    cdef int failed
    failed = cppm.readTriangleMesh(cmesh.meshlink.mesh,filebase,base);
    cppm.constructElementBoundaryElementsArray_triangle(cmesh.meshlink.mesh);
    failed = cppm.readTriangleElementBoundaryMaterialTypes(cmesh.meshlink.mesh,filebase,base);

cpdef void writeTriangleFiles(CMesh cmesh,
                             const char* filebase,
                             int base):
    cdef int failed
    failed = cppm.writeTriangleMesh(cmesh.meshlink.mesh,filebase,base);

cpdef void generateFromTetgenFiles(CMesh cmesh,
                                  const char* filebase,
                                  int base):
    cdef int failed
    failed = cppm.readTetgenMesh(cmesh.meshlink.mesh,filebase,base);
    cppm.constructElementBoundaryElementsArray_tetrahedron(cmesh.meshlink.mesh);
    failed = cppm.readTetgenElementBoundaryMaterialTypes(cmesh.meshlink.mesh,filebase,base);

cpdef void generateFromTetgenFilesParallel(CMesh cmesh,
                                          const char* filebase,
                                          int base):
    cdef int failed
    failed = cppm.readTetgenMesh(cmesh.meshlink.mesh,filebase,base);
    cppm.constructElementBoundaryElementsArray_tetrahedron(cmesh.meshlink.mesh);
    failed = cppm.readTetgenElementBoundaryMaterialTypes(cmesh.meshlink.mesh,filebase,base);

cpdef void writeTetgenFiles(CMesh cmesh,
                           const char* filebase,
                           int base):
    cdef int failed
    failed = cppm.writeTetgenMesh(cmesh.meshlink.mesh,filebase,base);

cpdef void write3dmFiles(CMesh cmesh,
                        const char* filebase,
                        int base):
    cdef int failed
    failed = cppm.write3dmMesh(cmesh.meshlink.mesh,filebase,base);

cpdef void write2dmFiles(CMesh cmesh,
                        const char* filebase,
                        int base):
    cdef int failed
    failed = cppm.write2dmMesh(cmesh.meshlink.mesh,filebase,base);

cpdef void generateFromHexFile(CMesh cmesh,
                              const char* filebase,
                              int base):
    cdef int failed
    failed = cppm.readHex(cmesh.meshlink.mesh,filebase,base);
    cppm.constructElementBoundaryElementsArray_hexahedron(cmesh.meshlink.mesh);

cpdef void generateFrom3DMFile(CMesh cmesh,
                              const char* filebase,
                              int base):
    cdef int failed
    failed = cppm.read3DM(cmesh.meshlink.mesh,filebase,base);
    cppm.constructElementBoundaryElementsArray_tetrahedron(cmesh.meshlink.mesh);
    failed = cppm.readBC(cmesh.meshlink.mesh,filebase,base);

cpdef void generateFrom2DMFile(CMesh cmesh,
                              const char* filebase,
                              int base):
    cdef int failed
    failed = cppm.read2DM(cmesh.meshlink.mesh,filebase,base);
    cppm.constructElementBoundaryElementsArray_triangle(cmesh.meshlink.mesh);
    failed = cppm.readBC(cmesh.meshlink.mesh,filebase,base);

cpdef void computeGeometricInfo_tetrahedron(CMesh cmesh):
    cppm.computeGeometricInfo_tetrahedron(cmesh.meshlink.mesh);

cpdef void allocateGeometricInfo_tetrahedron(CMesh cmesh):
    cppm.allocateGeometricInfo_tetrahedron(cmesh.meshlink.mesh);

cpdef void allocateNodeAndElementNodeDataStructures(CMesh cmesh,
                                                    int nElements_global,
                                                    int nNodes_global,
                                                    int nNodes_element):
    cppm.allocateNodeAndElementNodeDataStructures(cmesh.meshlink.mesh,nElements_global,nNodes_global,nNodes_element);

cpdef void constructElementBoundaryElementsArray(CMesh cmesh):
    if cmesh.meshlink.mesh.nNodes_element == 4:
        cppm.constructElementBoundaryElementsArray_tetrahedron(cmesh.meshlink.mesh);
    elif cmesh.meshlink.mesh.nNodes_element == 3:
        cppm.constructElementBoundaryElementsArray_triangle(cmesh.meshlink.mesh);
    else:
        cppm.constructElementBoundaryElementsArray_edge(cmesh.meshlink.mesh);

cpdef void generateTriangularMeshFromRectangularGrid(int nx,
                                                     int ny,
                                                     double Lx,
                                                     double Ly,
                                                     CMesh cmesh,
                                                     int triangleFlag):
    cppm.regularRectangularToTriangularMeshElements(nx,ny,cmesh.meshlink.mesh,triangleFlag);
    cppm.regularRectangularToTriangularMeshNodes(nx,ny,Lx,Ly,cmesh.meshlink.mesh);
    cppm.constructElementBoundaryElementsArray_triangle(cmesh.meshlink.mesh);
    cppm.regularRectangularToTriangularElementBoundaryMaterials(Lx,Ly,cmesh.meshlink.mesh);

cpdef void generateHexahedralMeshFromRectangularGrid(int nx,
                                                     int ny,
                                                     int nz,
                                                     int px,
                                                     int py,
                                                     int pz,
                                                     double Lx,
                                                     double Ly,
                                                     double Lz,
                                                     CMesh cmesh):
    cppm.regularHexahedralMeshElements(nx,ny,nz,px,py,pz,cmesh.meshlink.mesh);
    cppm.regularMeshNodes(nx,ny,nz,Lx,Ly,Lz,cmesh.meshlink.mesh);
    cppm.constructElementBoundaryElementsArray_hexahedron(cmesh.meshlink.mesh);
    cppm.regularHexahedralMeshElementBoundaryMaterials(Lx,Ly,Lz,cmesh.meshlink.mesh);

cpdef void generateQuadrilateralMeshFromRectangularGrid(int nx,
                                                        int ny,
                                                        int px,
                                                        int py,
                                                        double Lx,
                                                        double Ly,
                                                        CMesh cmesh):
    cppm.regularQuadrilateralMeshElements(nx,ny,cmesh.meshlink.mesh);
    cppm.regularMeshNodes2D(nx,ny,Lx,Ly,cmesh.meshlink.mesh);
    cppm.constructElementBoundaryElementsArray_quadrilateral(cmesh.meshlink.mesh);
    cppm.regularQuadrilateralMeshElementBoundaryMaterials(Lx,Ly,cmesh.meshlink.mesh);

cpdef void computeGeometricInfo_triangle(CMesh cmesh):
    cppm.computeGeometricInfo_triangle(cmesh.meshlink.mesh);

cpdef void allocateGeometricInfo_triangle(CMesh cmesh):
    cppm.allocateGeometricInfo_triangle(cmesh.meshlink.mesh);

cpdef void generateEdgeMeshFromRectangularGrid(int nx,
                                               double Lx,
                                               CMesh cmesh):
    cppm.edgeMeshElements(nx,cmesh.meshlink.mesh);
    cppm.regularEdgeMeshNodes(nx,Lx,cmesh.meshlink.mesh);
    cppm.constructElementBoundaryElementsArray_edge(cmesh.meshlink.mesh);

cpdef void computeGeometricInfo_edge(CMesh cmesh):
    cppm.computeGeometricInfo_edge(cmesh.meshlink.mesh);

cpdef void allocateGeometricInfo_edge(CMesh cmesh):
    cppm.allocateGeometricInfo_edge(cmesh.meshlink.mesh);


cpdef void computeGeometricInfo_hexahedron(CMesh cmesh):
    cppm.computeGeometricInfo_hexahedron(cmesh.meshlink.mesh);

cpdef void computeGeometricInfo_quadrilateral(CMesh cmesh):
    cppm.computeGeometricInfo_quadrilateral(cmesh.meshlink.mesh);

cpdef void allocateGeometricInfo_hexahedron(CMesh cmesh):
    cppm.allocateGeometricInfo_hexahedron(cmesh.meshlink.mesh);

cpdef void allocateGeometricInfo_quadrilateral(CMesh cmesh):
    cppm.allocateGeometricInfo_quadrilateral(cmesh.meshlink.mesh);

cpdef void computeGeometricInfo_NURBS(CMesh cmesh):
    cppm.computeGeometricInfo_NURBS(cmesh.meshlink.mesh);

cpdef void allocateGeometricInfo_NURBS(CMesh cmesh):
    cppm.allocateGeometricInfo_NURBS(cmesh.meshlink.mesh);
