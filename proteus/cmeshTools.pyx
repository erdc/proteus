# A type of -*- python -*- file
from cpython.cobject cimport PyCObject_AsVoidPtr
import cython
cimport cython
import numpy as np
cimport numpy as np
from libcpp cimport bool
from mpi4py import MPI
cimport mesh as cppm
from proteus import Comm

cdef class CMesh:
    def __init__(self):
        self.mesh = cppm.Mesh()
        cppm.initializeMesh(self.mesh)

    def buildPythonMeshInterface(self):
        cdef int dim1
        self.nElements_global = self.mesh.nElements_global
        self.nNodes_global = self.mesh.nNodes_global
        self.nNodes_element = self.mesh.nNodes_element
        self.nNodes_elementBoundary = self.mesh.nNodes_elementBoundary
        self.nElementBoundaries_element = self.mesh.nElementBoundaries_element
        self.nElementBoundaries_global = self.mesh.nElementBoundaries_global
        self.nInteriorElementBoundaries_global = self.mesh.nInteriorElementBoundaries_global
        self.nExteriorElementBoundaries_global = self.mesh.nExteriorElementBoundaries_global
        self.max_nElements_node = self.mesh.max_nElements_node
        self.nEdges_global = self.mesh.nEdges_global
        self.max_nNodeNeighbors_node = self.mesh.max_nNodeNeighbors_node
        self.elementNodesArray = np.asarray(<int[:self.mesh.nElements_global, :self.mesh.nNodes_element]> self.mesh.elementNodesArray)
        if self.mesh.nodeElementOffsets != NULL:
            self.nodeElementsArray = np.asarray(<int[:self.mesh.nodeElementOffsets[self.mesh.nNodes_global]]> self.mesh.nodeElementsArray)
            self.nodeElementOffsets = np.asarray(<int[:self.mesh.nNodes_global+1]> self.mesh.nodeElementOffsets)
        else:
            self.nodeElementsArray = np.empty(0, dtype=np.int32)
            self.nodeElementOffsets = np.empty(0, dtype=np.int32)

        self.elementNeighborsArray = np.asarray(<int[:self.mesh.nElements_global, :self.mesh.nElementBoundaries_element]> self.mesh.elementNeighborsArray)
        self.elementBoundariesArray = np.asarray(<int[:self.mesh.nElements_global, :self.mesh.nElementBoundaries_element]> self.mesh.elementBoundariesArray)
        self.elementBoundaryNodesArray = np.asarray(<int[:self.mesh.nElementBoundaries_global, :self.mesh.nNodes_elementBoundary]> self.mesh.elementBoundaryNodesArray)
        self.elementBoundaryElementsArray = np.asarray(<int[:self.mesh.nElementBoundaries_global, :2]> self.mesh.elementBoundaryElementsArray)
        self.elementBoundaryLocalElementBoundariesArray = np.asarray(<int[:self.mesh.nElementBoundaries_global, :2]> self.mesh.elementBoundaryLocalElementBoundariesArray)
        if self.mesh.nInteriorElementBoundaries_global:
            self.interiorElementBoundariesArray = np.asarray(<int[:self.mesh.nInteriorElementBoundaries_global]> self.mesh.interiorElementBoundariesArray)
        else:
            self.interiorElementBoundariesArray = np.empty(0, dtype=np.int32)
        self.exteriorElementBoundariesArray = np.asarray(<int[:self.mesh.nExteriorElementBoundaries_global]> self.mesh.exteriorElementBoundariesArray)
        self.edgeNodesArray = np.asarray(<int[:self.mesh.nEdges_global, :2]> self.mesh.edgeNodesArray)
        if self.mesh.nodeStarOffsets != NULL:
            self.nodeStarArray = np.asarray(<int[:self.mesh.nodeStarOffsets[self.mesh.nNodes_global]]> self.mesh.nodeStarArray)
        else:
            self.nodeStarArray = np.empty(0, dtype=np.int32)
        self.nodeStarOffsets = np.asarray(<int[:self.mesh.nNodes_global+1]> self.mesh.nodeStarOffsets)
        self.elementMaterialTypes = np.asarray(<int[:self.mesh.nElements_global]> self.mesh.elementMaterialTypes)
        self.elementBoundaryMaterialTypes = np.asarray(<int[:self.mesh.nElementBoundaries_global]> self.mesh.elementBoundaryMaterialTypes)
        self.nodeMaterialTypes = np.asarray(<int[:self.mesh.nNodes_global]> self.mesh.nodeMaterialTypes)
        self.nodeArray = np.asarray(<double[:self.mesh.nNodes_global, :3]> self.mesh.nodeArray)
        self.nx = self.mesh.nx
        self.ny = self.mesh.ny
        self.nz = self.mesh.nz
        self.px = self.mesh.px
        self.py = self.mesh.py
        self.pz = self.mesh.pz
        if self.mesh.elementIJK != NULL:
            self.elementIJK = np.asarray(<int[:self.mesh.nElements_global*3]> self.mesh.elementIJK)
        else:
            self.elementIJK = np.empty(0, dtype=np.int32)
        if self.mesh.weights != NULL:
            self.weights = np.asarray(<double[:self.mesh.nElements_global]> self.mesh.weights)
        else:
            self.weights = np.empty(0)
        if self.mesh.U_KNOT != NULL:
            self.U_KNOT = np.asarray(<double[:self.mesh.nx+self.mesh.px+1]> self.mesh.U_KNOT)
        else:
            self.U_KNOT = np.empty(0)
        if self.mesh.V_KNOT != NULL:
            self.V_KNOT = np.asarray(<double[:self.mesh.ny+self.mesh.py+1]> self.mesh.V_KNOT)
        else:
            self.V_KNOT = np.empty(0)
        if self.mesh.W_KNOT != NULL:
            self.W_KNOT = np.asarray(<double[:self.mesh.nz+self.mesh.pz+1]> self.mesh.W_KNOT)
        else:
            self.W_KNOT = np.empty(0)
        self.elementDiametersArray = np.asarray(<double[:self.mesh.nElements_global]> self.mesh.elementDiametersArray)
        self.elementInnerDiametersArray = np.asarray(<double[:self.mesh.nElements_global]> self.mesh.elementInnerDiametersArray)
        self.elementBoundaryDiametersArray = np.asarray(<double[:self.mesh.nElementBoundaries_global]> self.mesh.elementBoundaryDiametersArray)
        if self.mesh.elementBarycentersArray:
            self.elementBarycentersArray = np.asarray(<double[:self.mesh.nElements_global, :3]> self.mesh.elementBarycentersArray)
        if self.mesh.elementBoundaryBarycentersArray: 
            self.elementBoundaryBarycentersArray = np.asarray(<double[:self.mesh.nElementBoundaries_global, :3]> self.mesh.elementBoundaryBarycentersArray)
        if self.mesh.nodeDiametersArray:
            self.nodeDiametersArray = np.asarray(<double[:self.mesh.nNodes_global]> self.mesh.nodeDiametersArray)
        if self.mesh.nodeSupportArray:
            self.nodeSupportArray = np.asarray(<double[:self.mesh.nNodes_global]> self.mesh.nodeSupportArray)
        self.h = self.mesh.h
        self.hMin = self.mesh.hMin
        self.sigmaMax = self.mesh.sigmaMax
        self.volume = self.mesh.volume
    
    def buildPythonMeshInterfaceNoArrays(self):
        cdef int dim1
        self.nElements_global = self.mesh.nElements_global
        self.nNodes_global = self.mesh.nNodes_global
        self.nNodes_element = self.mesh.nNodes_element
        self.nNodes_elementBoundary = self.mesh.nNodes_elementBoundary
        self.nElementBoundaries_element = self.mesh.nElementBoundaries_element
        self.nElementBoundaries_global = self.mesh.nElementBoundaries_global
        self.nInteriorElementBoundaries_global = self.mesh.nInteriorElementBoundaries_global
        self.nExteriorElementBoundaries_global = self.mesh.nExteriorElementBoundaries_global
        self.max_nElements_node = self.mesh.max_nElements_node
        self.nEdges_global = self.mesh.nEdges_global
        self.max_nNodeNeighbors_node = self.mesh.max_nNodeNeighbors_node
        self.h = self.mesh.h
        self.hMin = self.mesh.hMin
        self.sigmaMax = self.mesh.sigmaMax
        self.volume = self.mesh.volume
    
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

def buildPythonMeshInterfaceNoArrays(cmesh):
    """
    function to be conform to old module, and to calls from MeshTools
    """
    cmesh.buildPythonMeshInterfaceNoArrays()
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
            cmesh.h,
            cmesh.hMin,
            cmesh.sigmaMax,
            cmesh.volume)

cdef CMesh CMesh_FromMesh(cppm.Mesh mesh):
    CMeshnew = CMesh()
    CMeshnew.mesh = mesh
    return CMeshnew

cdef class CMultilevelMesh:
    cdef cppm.MultilevelMesh multilevelMesh
    cdef public:
        int nLevels
        list cmeshList
        list elementParentsArrayList
        list elementChildrenArrayList
        list elementChildrenOffsetsList
    def __init__(self, CMesh cmesh, nLevels):
        self.multilevelMesh = cppm.MultilevelMesh()
        if cmesh.mesh.nNodes_element == 2:
            cppm.globallyRefineEdgeMesh(nLevels,
                                        cmesh.mesh,
                                        self.multilevelMesh,
                                        False);
            for i in range(1, nLevels):
                cppm.constructElementBoundaryElementsArray_edge(self.multilevelMesh.meshArray[i]);
                cppm.allocateGeometricInfo_edge(self.multilevelMesh.meshArray[i]);
                cppm.computeGeometricInfo_edge(self.multilevelMesh.meshArray[i]);
                cppm.assignElementBoundaryMaterialTypesFromParent(self.multilevelMesh.meshArray[i-1],
                                                                  self.multilevelMesh.meshArray[i],
                                                                  self.multilevelMesh.elementParentsArray[i],
                                                                  1);
        elif cmesh.mesh.nNodes_element == 3:
            cppm.globallyRefineTriangularMesh(nLevels,
                                              cmesh.mesh,
                                              self.multilevelMesh,
                                              False);
            for i in range(1, nLevels):
                cppm.constructElementBoundaryElementsArray_triangle(self.multilevelMesh.meshArray[i]);
                cppm.allocateGeometricInfo_triangle(self.multilevelMesh.meshArray[i]);
                cppm.computeGeometricInfo_triangle(self.multilevelMesh.meshArray[i]);
                cppm.assignElementBoundaryMaterialTypesFromParent(self.multilevelMesh.meshArray[i-1],
						                                     self.multilevelMesh.meshArray[i],
						                                     self.multilevelMesh.elementParentsArray[i],
						                                     2);
        elif cmesh.mesh.nNodes_element == 4 and cmesh.mesh.nNodes_elementBoundary == 2:
            cppm.globallyRefineQuadrilateralMesh(nLevels,
                                                 cmesh.mesh,
                                                 self.multilevelMesh,
                                                 False);
            for i in range(1, nLevels):
                cppm.constructElementBoundaryElementsArray_quadrilateral(self.multilevelMesh.meshArray[i]);
                cppm.allocateGeometricInfo_quadrilateral(self.multilevelMesh.meshArray[i]);
                cppm.computeGeometricInfo_quadrilateral(self.multilevelMesh.meshArray[i]);
                cppm.assignElementBoundaryMaterialTypesFromParent(self.multilevelMesh.meshArray[i-1],
						                                          self.multilevelMesh.meshArray[i],
						                                          self.multilevelMesh.elementParentsArray[i],
						                                          2);
        elif cmesh.mesh.nNodes_element == 4 and cmesh.mesh.nNodes_elementBoundary == 3:
            cppm.globallyRefineTetrahedralMesh(nLevels,
                                               cmesh.mesh,
                                               self.multilevelMesh,
                                               False);
            for i in range(1, nLevels):
                cppm.constructElementBoundaryElementsArray_tetrahedron(self.multilevelMesh.meshArray[i]);
                cppm.allocateGeometricInfo_tetrahedron(self.multilevelMesh.meshArray[i]);
                cppm.computeGeometricInfo_tetrahedron(self.multilevelMesh.meshArray[i]);
                cppm.assignElementBoundaryMaterialTypesFromParent(self.multilevelMesh.meshArray[i-1],
						                                     self.multilevelMesh.meshArray[i],
						                                     self.multilevelMesh.elementParentsArray[i],
						                                     3);
        elif cmesh.mesh.nNodes_element == 8 and cmesh.mesh.nNodes_elementBoundary == 4:
            cppm.globallyRefineHexahedralMesh(nLevels,
                                         cmesh.mesh,
                                         self.multilevelMesh,
                                         False);
            for i in range(1, nLevels):
                cppm.constructElementBoundaryElementsArray_hexahedron(self.multilevelMesh.
                                                                      meshArray[i]);
                cppm.allocateGeometricInfo_hexahedron(self.multilevelMesh.meshArray[i]);
                cppm.computeGeometricInfo_hexahedron(self.multilevelMesh.meshArray[i]);
                cppm.assignElementBoundaryMaterialTypesFromParent(self.multilevelMesh.meshArray[i-1],
						                                          self.multilevelMesh.meshArray[i],
						                                          self.multilevelMesh.elementParentsArray[i],
						                                          4);
        else:
            assert nLevels == 1, 'wrong nLevels'
            cppm.globallyRefineTetrahedralMesh(nLevels,
                                               cmesh.mesh,
                                               self.multilevelMesh,
                                               False);

    def buildPythonMultilevelMeshInterface(self):
        cdef int dim
        self.cmeshList = [CMesh_FromMesh(self.multilevelMesh.meshArray[0])]
        self.elementParentsArrayList = [np.zeros(0)]
        self.elementChildrenArrayList = []
        self.elementChildrenOffsetsList = []
        for n in range(1, self.multilevelMesh.nLevels):
            self.cmeshList += [CMesh_FromMesh(self.multilevelMesh.meshArray[n])]
            dim = self.multilevelMesh.meshArray[n].nElements_global
            self.elementParentsArrayList += [np.asarray(<int[:dim]> self.multilevelMesh.elementParentsArray[n])]
            dim = self.multilevelMesh.elementChildrenOffsets[n-1][self.multilevelMesh.meshArray[n-1].nElements_global]
            self.elementChildrenArrayList += [np.asarray(<int[:dim]> self.multilevelMesh.elementChildrenArray[n-1])]
            dim = self.multilevelMesh.meshArray[n-1].nElements_global+1;
            self.elementChildrenOffsetsList += [np.asarray(<int[:dim]> self.multilevelMesh.elementChildrenOffsets[n-1])]

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
    cppm.regularHexahedralToTetrahedralMeshElements(nx,ny,nz,cmesh.mesh);
    cppm.regularHexahedralToTetrahedralMeshNodes(nx,ny,nz,Lx,Ly,Lz,cmesh.mesh);
    cppm.constructElementBoundaryElementsArray_tetrahedron(cmesh.mesh);
    cppm.regularHexahedralToTetrahedralElementBoundaryMaterials(Lx,Ly,Lz,cmesh.mesh);

cpdef void cmeshToolsComputeGeometricInfo_tetrahedron(CMesh cmesh):
    cppm.computeGeometricInfo_tetrahedron(cmesh.mesh)

# cdef static void cmeshToolsLocallyRefineMultilevelMesh(CMultilevelMesh cmesh,


def generateFromTriangleMesh(CMesh cmesh,
                             object trimesh,
			     int base):
    cppm.setFromTriangleElements(<cppm.triangulateio*>(PyCObject_AsVoidPtr(trimesh)),
                                 cmesh.mesh,
                                 base)
    cppm.setFromTriangleNodes(<cppm.triangulateio*>(PyCObject_AsVoidPtr(trimesh)),
                              cmesh.mesh,
                              base)
    cppm.constructElementBoundaryElementsArray_triangle(cmesh.mesh)
    cppm.copyElementBoundaryMaterialTypesFromTriangle(<cppm.triangulateio*>(PyCObject_AsVoidPtr(trimesh)),
					              cmesh.mesh,
                                                      base)

cpdef void generateFromTriangleFiles(CMesh cmesh,
                                    const char* filebase,
                                    int base):

    cdef int failed
    failed = cppm.readTriangleMesh(cmesh.mesh,filebase,base);
    cppm.constructElementBoundaryElementsArray_triangle(cmesh.mesh);
    failed = cppm.readTriangleElementBoundaryMaterialTypes(cmesh.mesh,filebase,base);

cpdef void writeTriangleFiles(CMesh cmesh,
                             const char* filebase,
                             int base):
    cdef int failed
    failed = cppm.writeTriangleMesh(cmesh.mesh,filebase,base);

cpdef void generateFromTetgenFiles(CMesh cmesh,
                                  const char* filebase,
                                  int base):
    cdef int failed
    failed = cppm.readTetgenMesh(cmesh.mesh,filebase,base);
    cppm.constructElementBoundaryElementsArray_tetrahedron(cmesh.mesh);
    failed = cppm.readTetgenElementBoundaryMaterialTypes(cmesh.mesh,filebase,base);

cpdef void generateFromTetgenFilesParallel(CMesh cmesh,
                                          const char* filebase,
                                          int base):
    cdef int failed
    failed = cppm.readTetgenMesh(cmesh.mesh,filebase,base);
    cppm.constructElementBoundaryElementsArray_tetrahedron(cmesh.mesh);
    failed = cppm.readTetgenElementBoundaryMaterialTypes(cmesh.mesh,filebase,base);

cpdef void writeTetgenFiles(CMesh cmesh,
                           const char* filebase,
                           int base):
    cdef int failed
    failed = cppm.writeTetgenMesh(cmesh.mesh,filebase,base);

cpdef void write3dmFiles(CMesh cmesh,
                        const char* filebase,
                        int base):
    cdef int failed
    failed = cppm.write3dmMesh(cmesh.mesh,filebase,base);

cpdef void write2dmFiles(CMesh cmesh,
                        const char* filebase,
                        int base):
    cdef int failed
    failed = cppm.write2dmMesh(cmesh.mesh,filebase,base);

cpdef void generateFromHexFile(CMesh cmesh,
                              const char* filebase,
                              int base):
    cdef int failed
    failed = cppm.readHex(cmesh.mesh,filebase,base);
    cppm.constructElementBoundaryElementsArray_hexahedron(cmesh.mesh);

cpdef void generateFrom3DMFile(CMesh cmesh,
                              const char* filebase,
                              int base):
    cdef int failed
    failed = cppm.read3DM(cmesh.mesh,filebase,base);
    cppm.constructElementBoundaryElementsArray_tetrahedron(cmesh.mesh);
    failed = cppm.readBC(cmesh.mesh,filebase,base);

cpdef void generateFrom2DMFile(CMesh cmesh,
                              const char* filebase,
                              int base):
    cdef int failed
    failed = cppm.read2DM(cmesh.mesh,filebase,base);
    cppm.constructElementBoundaryElementsArray_triangle(cmesh.mesh);
    failed = cppm.readBC(cmesh.mesh,filebase,base);

cpdef void computeGeometricInfo_tetrahedron(CMesh cmesh):
    cppm.computeGeometricInfo_tetrahedron(cmesh.mesh);

cpdef void allocateGeometricInfo_tetrahedron(CMesh cmesh):
    cppm.allocateGeometricInfo_tetrahedron(cmesh.mesh);

cpdef void allocateNodeAndElementNodeDataStructures(CMesh cmesh,
                                                    int nElements_global,
                                                    int nNodes_global,
                                                    int nNodes_element):
    cppm.allocateNodeAndElementNodeDataStructures(cmesh.mesh,nElements_global,nNodes_global,nNodes_element);

cpdef void constructElementBoundaryElementsArray(CMesh cmesh):
    if cmesh.mesh.nNodes_element == 4:
        cppm.constructElementBoundaryElementsArray_tetrahedron(cmesh.mesh);
    elif cmesh.mesh.nNodes_element == 3:
        cppm.constructElementBoundaryElementsArray_triangle(cmesh.mesh);
    else:
        cppm.constructElementBoundaryElementsArray_edge(cmesh.mesh);

cpdef void generateTriangularMeshFromRectangularGrid(int nx,
                                                     int ny,
                                                     double Lx,
                                                     double Ly,
                                                     CMesh cmesh,
                                                     int triangleFlag):
    cppm.regularRectangularToTriangularMeshElements(nx,ny,cmesh.mesh,triangleFlag);
    cppm.regularRectangularToTriangularMeshNodes(nx,ny,Lx,Ly,cmesh.mesh);
    cppm.constructElementBoundaryElementsArray_triangle(cmesh.mesh);
    cppm.regularRectangularToTriangularElementBoundaryMaterials(Lx,Ly,cmesh.mesh);

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
    cppm.regularHexahedralMeshElements(nx,ny,nz,px,py,pz,cmesh.mesh);
    cppm.regularMeshNodes(nx,ny,nz,Lx,Ly,Lz,cmesh.mesh);
    cppm.constructElementBoundaryElementsArray_hexahedron(cmesh.mesh);
    cppm.regularHexahedralMeshElementBoundaryMaterials(Lx,Ly,Lz,cmesh.mesh);

cpdef void generateQuadrilateralMeshFromRectangularGrid(int nx,
                                                        int ny,
                                                        int px,
                                                        int py,
                                                        double Lx,
                                                        double Ly,
                                                        CMesh cmesh):
    cppm.regularQuadrilateralMeshElements(nx,ny,cmesh.mesh);
    cppm.regularMeshNodes2D(nx,ny,Lx,Ly,cmesh.mesh);
    cppm.constructElementBoundaryElementsArray_quadrilateral(cmesh.mesh);
    cppm.regularQuadrilateralMeshElementBoundaryMaterials(Lx,Ly,cmesh.mesh);

cpdef void computeGeometricInfo_triangle(CMesh cmesh):
    cppm.computeGeometricInfo_triangle(cmesh.mesh);

cpdef void allocateGeometricInfo_triangle(CMesh cmesh):
    cppm.allocateGeometricInfo_triangle(cmesh.mesh);

cpdef void generateEdgeMeshFromRectangularGrid(int nx,
                                               double Lx,
                                               CMesh cmesh):
    cppm.edgeMeshElements(nx,cmesh.mesh);
    cppm.regularEdgeMeshNodes(nx,Lx,cmesh.mesh);
    cppm.constructElementBoundaryElementsArray_edge(cmesh.mesh);

cpdef void computeGeometricInfo_edge(CMesh cmesh):
    cppm.computeGeometricInfo_edge(cmesh.mesh);

cpdef void allocateGeometricInfo_edge(CMesh cmesh):
    cppm.allocateGeometricInfo_edge(cmesh.mesh);


cpdef void computeGeometricInfo_hexahedron(CMesh cmesh):
    cppm.computeGeometricInfo_hexahedron(cmesh.mesh);

cpdef void computeGeometricInfo_quadrilateral(CMesh cmesh):
    cppm.computeGeometricInfo_quadrilateral(cmesh.mesh);

cpdef void allocateGeometricInfo_hexahedron(CMesh cmesh):
    cppm.allocateGeometricInfo_hexahedron(cmesh.mesh);

cpdef void allocateGeometricInfo_quadrilateral(CMesh cmesh):
    cppm.allocateGeometricInfo_quadrilateral(cmesh.mesh);

cpdef void computeGeometricInfo_NURBS(CMesh cmesh):
    cppm.computeGeometricInfo_NURBS(cmesh.mesh);

cpdef void allocateGeometricInfo_NURBS(CMesh cmesh):
    cppm.allocateGeometricInfo_NURBS(cmesh.mesh);

def locallyRefineMultilevelMesh(int nSpace,
                                          CMultilevelMesh cmultilevelMesh,
                                          np.ndarray elementTagArray,
                                          int refineTypeFlag=0):
    cdef int failed,finestLevel
    if nSpace == 1:
        failed = cppm.locallyRefineEdgeMesh(cmultilevelMesh.multilevelMesh,
                                            <int*>(elementTagArray.data))
        finestLevel = cmultilevelMesh.multilevelMesh.nLevels
        cppm.constructElementBoundaryElementsArray_edge(cmultilevelMesh.multilevelMesh.meshArray[finestLevel-1])
        cppm.allocateGeometricInfo_edge(cmultilevelMesh.multilevelMesh.meshArray[finestLevel-1])
        cppm.computeGeometricInfo_edge(cmultilevelMesh.multilevelMesh.meshArray[finestLevel-1])
        if finestLevel > 1:
            cppm.assignElementBoundaryMaterialTypesFromParent(cmultilevelMesh.multilevelMesh.meshArray[finestLevel-2],
                                                              cmultilevelMesh.multilevelMesh.meshArray[finestLevel-1],
                                                              cmultilevelMesh.multilevelMesh.elementParentsArray[finestLevel-1],
                                                              1)
    elif nSpace == 2:
        if refineTypeFlag==1:
            failed = cppm.locallyRefineTriangleMesh_4T(cmultilevelMesh.multilevelMesh,
                                                       <int*>(elementTagArray.data))
        elif refineTypeFlag==2:
            failed = cppm.locallyRefineTriangleMesh_redGreen(cmultilevelMesh.multilevelMesh,<int*>(elementTagArray.data))
        else:
            failed = cppm.locallyRefineTriangleMesh(cmultilevelMesh.multilevelMesh,
                                                    <int*>(elementTagArray.data))
        finestLevel = cmultilevelMesh.multilevelMesh.nLevels
        cppm.constructElementBoundaryElementsArray_triangle(cmultilevelMesh.multilevelMesh.meshArray[finestLevel-1])
        cppm.allocateGeometricInfo_triangle(cmultilevelMesh.multilevelMesh.meshArray[finestLevel-1])
        cppm.computeGeometricInfo_triangle(cmultilevelMesh.multilevelMesh.meshArray[finestLevel-1])
        if finestLevel > 1:
            cppm.assignElementBoundaryMaterialTypesFromParent(cmultilevelMesh.multilevelMesh.meshArray[finestLevel-2],
                                                              cmultilevelMesh.multilevelMesh.meshArray[finestLevel-1],
                                                              cmultilevelMesh.multilevelMesh.elementParentsArray[finestLevel-1],
                                                              2)
    else:
        print("locallyRefine nSpace= {0:d} not implemented! Returning.".format(nSpace))

def setNewestNodeBases(int nSpace, CMultilevelMesh cmultilevelMesh):
    if nSpace == 2:
        failed = cppm.setNewestNodeBasesToLongestEdge(cmultilevelMesh.multilevelMesh)
    else:
        print("setNewestNodeBases= {0:d} not implemented! Returning".format(nSpace))
