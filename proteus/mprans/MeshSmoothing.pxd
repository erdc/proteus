cimport cython
import numpy as np
cimport numpy as np
from libcpp cimport bool

cdef void cySmoothNodesLaplace(double[:,:] nodeArray_,
                               int[:] nodeStarOffsets,
                               int[:] nodeStarArray,
                               int[:] nodeMaterialTypes,
                               int nNodes_owned,
                               int nd,
                               bool simultaneous=*,
                               bool smoothBoundaries=*,
                               int[:] fixedNodesBoolArray=*,
                               double alpha=*)

cdef void cySmoothNodesCentroid(double[:,:] nodeArray_,
                                int[:] nodeElementOffsets,
                                int[:] nodeElementsArray,
                                int[:] nodeMaterialTypes,
                                double[:] elementVolumesArray,
                                double[:,:] elementBarycentersArray,
                                int[:,:] elementNodesArray,
                                int nNodes_owned,
                                int[:] fixedNodesBoolArray,
                                bool simultaneous=*,
                                bool smoothBoundaries=*,
                                double alpha=*)

cdef void cyUpdateDilationElements(double[:] elementDilationArray_,
                                   double[:] elementVolumeArray,
                                   double[:] elementVolumeTargetArray,
                                   int nElements)

cdef void cyUpdateDistortionElements(double[:] elementDistortionArray_,
                                     double[:,:,:,:] J_array,
                                     double[:,:] detJ_array,
                                     int nd,
                                     int nElements)

cdef void cyUpdateInverseMeanRatioTriangleNodes(double[:] IMRNodesArray_,
                                                double[:,:] nodeArray,
                                                int[:,:] elementNodesArray,
                                                int[:] nodeElementOffsets,
                                                int[:] nodeElementsArray,
                                                int nElements,
                                                int nNodes,
                                                bool el_average=*)

cdef void cyUpdateInverseMeanRatioTriangleElements(double[:] IMRElementsArray_,
                                                   double[:,:] nodeArray,
                                                   int[:,:] elementNodesArray,
                                                   int nElements)

cdef double cyGetInverseMeanRatioSingleTriangle(int node0,
                                                double[:,:] nodeArray,
                                                int[:,:] elementNodesArray,
                                                int[:] nodeElementOffsets,
                                                int[:] nodeElementsArray,
                                                bool el_average=*)

cdef double[:,:] cySmoothNodesQuality(double[:] distortion,
                                      double[:] dilation,
                                      double[:,:] nodeArray,
                                      int nNodes_owned,
                                      int[:] nodeMaterialTypes,
                                      int[:] nodeElementOffsets,
                                      int[:] nodeElementsArray,
                                      int[:,:] elementNodesArray,
                                      bool apply_directly=*)

cdef int pyxGetLocalNearestNode(double[:] coords,
                                double[:,:] nodeArray,
                                int[:] nodeStarOffsets,
                                int[:] nodeStarArray,
                                int node)

cdef int pyxGetLocalNearestElement(double[:] coords,
                                   double[:,:] elementBarycentersArray,
                                   int[:,:] elementNeighborsArray,
                                   int eN)

cdef int[:] pyxGetLocalNearestElementIntersection(double[:] coords,
                                                  double[:] starting_coords,
                                                  double[:,:,:] elementBoundaryNormalsArray,
                                                  int[:,:] elementBoundariesArray,
                                                  double[:,:] elementBoundaryBarycentersArray,
                                                  int[:,:] elementBoundaryElementsArray,
                                                  int[:] exteriorElementBoundariesBoolArray,
                                                  int eN)

cdef int pyxGetLocalNearestElementAroundNode(double[:] coords,
                                             int[:] nodeElementOffsets,
                                             int[:] nodeElementsArray,
                                             double[:,:] elementBarycentersArray,
                                             int node)

cdef void pyxUpdateElementBoundaryNormalsTetra(double[:,:,:] elementBoundaryNormalsArray_,
                                               double[:,:] nodeArray,
                                               int[:,:] elementBoundariesArray,
                                               int[:,:] elementBoundaryNodesArray,
                                               double[:,:] elementBoundaryBarycentersArray,
                                               double[:,:] elementBarycentersArray,
                                               int nElements)

cdef void pyxUpdateElementBoundaryNormalsTriangle(double[:,:,:] elementBoundaryNormalsArray_,
                                                  double[:,:] nodeArray,
                                                  int[:,:] elementBoundariesArray,
                                                  int[:,:] elementBoundaryNodesArray,
                                                  double[:,:] elementBoundaryBarycentersArray,
                                                  double[:,:] elementBarycentersArray,
                                                  int nElements)

cdef void cyUpdateElementVolumesTriangle(double[:] elementVolumesArray_,
                                         int[:,:] elementNodesArray,
                                         double[:,:] nodeArray,
                                         int nElements)

cdef void cyUpdateElementVolumesTetra(double[:] elementVolumesArray_,
                                      int[:,:] elementNodesArray,
                                      double[:,:] nodeArray,
                                      int nElements)

cdef void cyUpdateElementBarycenters(double[:,:] elementBarycentersArray_,
                                     int[:,:] elementNodesArray,
                                     double[:,:] nodeArray,
                                     int nElements)

cdef np.ndarray cyGetCornerNodesTriangle(double[:,:] nodeArray,
                                         int[:] nodeStarArray,
                                         int[:] nodeStarOffsets,
                                         int[:] nodeMaterialTypes,
                                         int nNodes)

cdef int[:] cyCheckOwnedVariable(int variable_nb_local,
                                 int rank,
                                 int nVariables_owned,
                                 int[:] variableNumbering_subdomain2global,
                                 int[:] variableOffsets_subdomain_owned)

cdef int[:] cyGetGlobalVariable(int variable_nb_local,
                                int nVariables_owned,
                                int[:] variableNumbering_subdomain2global,
                                int[:] variableOffsets_subdomain_owned)

cdef int cyGetLocalVariable(int variable_nb_global,
                            int rank,
                            int nVariables_owned,
                            int[:] variableNumbering_subdomain2global,
                            int[:] variableOffsets_subdomain_owned)

cdef double[:] cyScalarRecoveryAtNodes(double[:] scalars,
                                       int[:] nodeElementsArray,
                                       int[:] nodeElementOffsets)

cdef double[:] cyScalarRecoveryAtNodesWeighted(double[:] scalars,
                                               int[:] nodeElementsArray,
                                               int[:] nodeElementOffsets,
                                               double[:] detJ_array,
                                               int nNodes)

cdef double[:,:] cyVectorRecoveryAtNodesWeighted(double[:,:] vectors,
                                                 int[:] nodeElementsArray,
                                                 int[:] nodeElementOffsets,
                                                 double[:,:] detJ_array,
                                                 int nd)

cdef double[:,:] cyVectorRecoveryAtNodes(double[:,:] vectors,
                                         int[:] nodeElementsArray,
                                         int[:] nodeElementOffsets,
                                         int nd)
