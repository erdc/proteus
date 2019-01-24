from libcpp cimport bool

cdef extern from "mesh.h":
    cdef struct Mesh:
      int nElements_global
      int nNodes_global
      int nNodes_element
      int nNodes_elementBoundary
      int nElementBoundaries_element
      int nElementBoundaries_global
      int nInteriorElementBoundaries_global
      int nExteriorElementBoundaries_global
      int max_nElements_node
      int nEdges_global
      int max_nNodeNeighbors_node

      int *elementNodesArray
      int *nodeElementsArray
      int *nodeElementOffsets
      int *elementNeighborsArray
      int *elementBoundariesArray
      int *elementBoundaryNodesArray
      int *elementBoundaryElementsArray
      int *elementBoundaryLocalElementBoundariesArray
      int *interiorElementBoundariesArray
      int *exteriorElementBoundariesArray
      int *edgeNodesArray
      int *nodeStarArray
      int *nodeStarOffsets
      int *elementMaterialTypes
      int *elementBoundaryMaterialTypes
      int *nodeMaterialTypes

      int *elementIJK
      double *weights
      double *U_KNOT,*V_KNOT,*W_KNOT
      int nx,ny,nz
      int px,py,pz

      double *nodeArray,*elementDiametersArray,*elementInnerDiametersArray,*elementBoundaryDiametersArray
      double *elementBarycentersArray, *elementBoundaryBarycentersArray
      double *nodeDiametersArray,*nodeSupportArray
      double h,hMin,sigmaMax,volume
      int * newestNodeBases

      int *elementOffsets_subdomain_owned,
      int *elementNumbering_subdomain2global,
      int *nodeOffsets_subdomain_owned,
      int *nodeNumbering_subdomain2global,
      int *elementBoundaryOffsets_subdomain_owned,
      int *elementBoundaryNumbering_subdomain2global,
      int *edgeOffsets_subdomain_owned,
      int *edgeNumbering_subdomain2global
      Mesh* subdomainp

    cdef inline void initializeMesh(Mesh& mesh)

    cdef inline void deleteMesh(Mesh& mesh)

    cdef struct MultilevelMesh:
        int nLevels
        Mesh* meshArray
        int **elementParentsArray
        int **elementChildrenArray
        int **elementChildrenOffsets

    cdef inline void initializeMultilevelMesh(MultilevelMesh& multilevelMesh)

    cdef inline void deleteMultilevelMesh(MultilevelMesh& multilevelMesh)

    cdef int edgeMeshElements(const int& nx,
                              Mesh& mesh)
    cdef int regularEdgeMeshNodes(const int& nx,
                                  const double& Lx,
                                  Mesh& mesh)
    cdef int globallyRefineEdgeMesh(const int& nLevels,
                                    Mesh& mesh,
                                    MultilevelMesh& multilevelMesh)
    cdef int globallyRefineEdgeMesh(const int& nLevels,
                                    Mesh& mesh,
                                    MultilevelMesh& multilevelMesh,
                                    bool averageNewNodeFlags)
    cdef int locallyRefineEdgeMesh(MultilevelMesh& multilevelMesh,
			                       int * elementTagArray)
    cdef int locallyRefineTriangleMesh(MultilevelMesh& multilevelMesh,
				                       int * elementTagArray)
    cdef int locallyRefineTriangleMesh_4T(MultilevelMesh& multilevelMesh,
				                          int * elementTagArray)
    cdef int locallyRefineTriangleMesh_redGreen(MultilevelMesh& multilevelMesh,
				                              int * elementTagArray)
    cdef int setNewestNodeBasesToLongestEdge(MultilevelMesh& multilevelMesh)

    cdef int regularRectangularToTriangularMeshElements(const int& nx,
                                                        const int& ny,
                                                        Mesh& mesh,
                                                        int triangleFlag)
    cdef int regularRectangularToTriangularMeshNodes(const int& nx,
                                                     const int& ny,
                                                     const double& Lx,
                                                     const double& Ly,
                                                     Mesh& mesh)
    cdef int regularRectangularToTriangularElementBoundaryMaterials(const double& Lx,
                                                                    const double& Ly,
                                                                    Mesh& mesh)
    cdef int globallyRefineTriangularMesh(const int& nLevels,
                                          Mesh& mesh,
                                          MultilevelMesh& multilevelMesh)
    cdef int globallyRefineTriangularMesh(const int& nLevels,
                                          Mesh& mesh,
                                          MultilevelMesh& multilevelMesh,
                                          bool averageNewNodeFlags)

    cdef int regularQuadrilateralMeshElements(const int& nx,
                                              const int& ny,
                                              Mesh& mesh)
    cdef int regularQuadrilateralMeshElementBoundaryMaterials(const double& Lx,
                                                              const double& Ly,
                                                              Mesh& mesh)
    cdef int globallyRefineQuadrilateralMesh(const int& nLevels,
                                             Mesh& mesh,
                                             MultilevelMesh& multilevelMesh)
    cdef int globallyRefineQuadrilateralMesh(const int& nLevels,
                                             Mesh& mesh,
                                             MultilevelMesh& multilevelMesh,
                                             bool averageNewNodeFlags)

    cdef int regularMeshNodes(const int& nx,
                              const int& ny,
                              const int& nz,
                              const double& Lx,
                              const double& Ly,
                              const double& Lz,
                              Mesh& mesh)
    cdef int regularMeshNodes2D(const int& nx,
                                const int& ny,
                                const double& Lx,
                                const double& Ly,
                                Mesh& mesh)
    cdef int regularHexahedralMeshElementBoundaryMaterials(const double& Lx,
                                                           const double& Ly,
                                                           const double& Lz,
                                                           Mesh& mesh)
    cdef int regularHexahedralToTetrahedralMeshNodes(const int& nx,
                                                     const int& ny,
                                                     const int& nz,
                                                     const double& Lx,
                                                     const double& Ly,
                                                     const double& Lz,
                                                     Mesh& mesh)
    cdef int regularHexahedralToTetrahedralMeshElements(const int& nx,
                                                        const int& ny,
                                                        const int& nz,
                                                        Mesh& mesh)
    cdef int regularHexahedralToTetrahedralElementBoundaryMaterials(const double& Lx,
                                                                    const double& Ly,
                                                                    const double& Lz,
                                                                    Mesh& mesh)
    cdef int regularHexahedralMeshElements(const int& nx,
                                           const int& ny,
                                           const int& nz,
                                           const int& px,
                                           const int& py,
                                           const int& pz,
                                           Mesh& mesh)
    cdef int regularNURBSMeshElements(const int& nx,
                                      const int& ny,
                                      const int& nz,
                                      const int& px,
                                      const int& py,
                                      const int& pz,
                                      Mesh& mesh)
    cdef int globallyRefineHexahedralMesh(const int& nLevels,
                                          Mesh& mesh,
                                          MultilevelMesh& multilevelMesh)
    cdef int globallyRefineHexahedralMesh(const int& nLevels,
                                          Mesh& mesh,
                                          MultilevelMesh& multilevelMesh,
                                          bool averageNewNodeFlags)

    cdef int globallyRefineTetrahedralMesh(const int& nLevels,
                                           Mesh& mesh,
                                           MultilevelMesh& multilevelMesh)
    cdef int globallyRefineTetrahedralMesh(const int& nLevels,
                                           Mesh& mesh,
                                           MultilevelMesh& multilevelMesh,
                                           bool averageNewNodeFlags)

    cdef int constructElementBoundaryElementsArray_edge(Mesh& mesh)
    cdef int constructElementBoundaryElementsArray_triangle(Mesh& mesh)
    cdef int constructElementBoundaryElementsArray_quadrilateral(Mesh& mesh)
    cdef int constructElementBoundaryElementsArray_tetrahedron(Mesh& mesh)
    cdef int constructElementBoundaryElementsArray_hexahedron(Mesh& mesh)
    cdef int constructElementBoundaryElementsArray_NURBS(Mesh& mesh)

    cdef int constructElementBoundaryElementsArrayWithGivenElementBoundaryNumbers_edge(Mesh& mesh)
    cdef int constructElementBoundaryElementsArrayWithGivenElementBoundaryNumbers_triangle(Mesh& mesh)
    cdef int constructElementBoundaryElementsArrayWithGivenElementBoundaryNumbers_quadrilateral(Mesh& mesh)
    cdef int constructElementBoundaryElementsArrayWithGivenElementBoundaryNumbers_tetrahedron(Mesh& mesh)

    cdef int constructElementBoundaryElementsArrayWithGivenElementBoundaryAndEdgeNumbers_edge(Mesh& mesh)
    cdef int constructElementBoundaryElementsArrayWithGivenElementBoundaryAndEdgeNumbers_triangle(Mesh& mesh)
    cdef int constructElementBoundaryElementsArrayWithGivenElementBoundaryAndEdgeNumbers_quadrilateral(Mesh& mesh)
    cdef int constructElementBoundaryElementsArrayWithGivenElementBoundaryAndEdgeNumbers_tetrahedron(Mesh& mesh)
    cdef int constructElementBoundaryElementsArrayWithGivenElementBoundaryAndEdgeNumbers_hexahedron(Mesh& mesh)
    cdef int constructElementBoundaryElementsArrayWithGivenElementBoundaryAndEdgeNumbers_NURBS(Mesh& mesh)

    # cdef int writeElements(std::ostream& meshFile,
    #                        const Mesh& mesh)
    # cdef int writeNodes(std::ostream& meshFile,
    #                     const Mesh& mesh)
    # cdef int readElements(std::istream& meshFile,
    #                       Mesh& mesh)


    cdef int allocateGeometricInfo_tetrahedron(Mesh& mesh)
    cdef int allocateGeometricInfo_triangle(Mesh& mesh)
    cdef int allocateGeometricInfo_edge(Mesh& mesh)
    cdef int allocateGeometricInfo_quadrilateral(Mesh& mesh)
    cdef int allocateGeometricInfo_hexahedron(Mesh& mesh)
    cdef int allocateGeometricInfo_NURBS(Mesh& mesh)

    cdef int computeGeometricInfo_tetrahedron(Mesh& mesh)
    cdef int computeGeometricInfo_triangle(Mesh& mesh)
    cdef int computeGeometricInfo_edge(Mesh& mesh)
    cdef int computeGeometricInfo_hexahedron(Mesh& mesh)
    cdef int computeGeometricInfo_quadrilateral(Mesh& mesh)
    cdef int computeGeometricInfo_NURBS(Mesh& mesh)

    cdef int assignElementBoundaryMaterialTypesFromParent(Mesh& parentMesh,
                                                          Mesh& childMesh,
                                                          const int* levelElementParentsArray,
						                                  const int& nSpace_global)
    cdef int allocateNodeAndElementNodeDataStructures(Mesh& mesh,
                                                      int nElements_global,
                                                      int nNodes_global,
                                                      int nNodes_element)
    cdef struct triangulateio

    cdef int setFromTriangleElements(triangulateio* trimesh,
                                     Mesh& mesh,
                                     int base)
    cdef int setFromTriangleNodes(triangulateio* trimesh,
                                  Mesh& mesh,
                                  int base)
    cdef int readTriangleMesh(Mesh& mesh,
                              const char* filebase,
                              int base)
    cdef int readTriangleElementBoundaryMaterialTypes(Mesh& mesh,
                                                      const char* filebase,
                                                      int base)
    cdef int writeTriangleMesh(Mesh& mesh,
                               const char* filebase,
                               int base)
    cdef int readTetgenMesh(Mesh& mesh,
                            const char* filebase,
                            int base)
    cdef int readTetgenElementBoundaryMaterialTypes(Mesh& mesh,
                                                    const char* filebase,
                                                    int base)
    cdef int writeTetgenMesh(Mesh& mesh,
                             const char* filebase,
                             int base)
    cdef int read3DM(Mesh& mesh,
                     const char* filebase,
                     int indexBase)
    cdef int read2DM(Mesh& mesh,
                     const char* filebase,
                     int indexBase)
    cdef int readHex(Mesh& mesh,
                     const char* filebase,
                     int indexBase)
    cdef int readBC(Mesh& mesh,
                    const char* filebase,
                    int indexBase)
    cdef int write3dmMesh(Mesh& mesh,
                          const char* filebase,
                          int base)
    cdef int write2dmMesh(Mesh& mesh,
                          const char* filebase,
                          int base)
    cdef int copyElementBoundaryMaterialTypesFromTriangle(triangulateio* trimesh,
						                                  Mesh& mesh,
                                                          int base)
