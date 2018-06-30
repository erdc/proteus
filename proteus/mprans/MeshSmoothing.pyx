import numpy as np
cimport numpy as np
from libcpp cimport bool

def smoothNodesLaplace(double[:,:] nodeArray_,
                       int[:] nodeStarOffsets,
                       int[:] nodeStarArray,
                       int[:] nodeMaterialTypes,
                       int nNodes_owned,
                       bool simultaneous=False,
                       bool smoothBoundaries=True,
                       int[:] fixedNodesBoolArray=None,
                       double alpha=0.):
    """
    Laplace Smoothing:
    Mesh nodes are displaced to the centroid of neighbouring nodes

    Parameters
    ----------
    nodeArray_: double[:,:]
        array of mesh nodes
        (!) will be modified with smoothing
    nodeStarOffsets: int[:]
        array to get offsets for neighbouring node numbers
    nodeStarArray: int[:]
        array of neighbouring nodes
    nodeMaterialTypes: int[:]
        array to know which node is interior node (materialType=0)
    nNodes_owned: int
        number of nodes owned
    simultaneous: bool
        if True: simultaneous smoothing
        if False: sequential smoothing
    smoothBoundaries: bool
        if True: boundaries are smoothed (only with surrounding boundary nodes)
        if False: boudnary nodes are fixed
    fixedNodesBoolArray: int[:]
        array of same length as nodeArray with:
        fixedNodesBoolArray[node] = 1 -> node is fixed
        fixedNodesBoolArray[node] = 0 -> node is not fixed
    alpha: double
        penalty term related to original position of nodes (0 <= alpha <= 1)
    """
    cySmoothNodesLaplace(nodeArray_=nodeArray_,
                         nodeStarOffsets=nodeStarOffsets,
                         nodeStarArray=nodeStarArray,
                         nodeMaterialTypes=nodeMaterialTypes,
                         nNodes_owned=nNodes_owned,
                         simultaneous=simultaneous,
                         smoothBoundaries=smoothBoundaries,
                         fixedNodesBoolArray=fixedNodesBoolArray,
                         alpha=alpha)

def smoothNodesCentroid(double[:,:] nodeArray_,
                        int[:] nodeElementOffsets,
                        int[:] nodeElementsArray,
                        int[:] nodeMaterialTypes,
                        double[:,:] elementBarycentersArray,
                        double[:] elementVolumesArray,
                        int[:,:] elementNodesArray,
                        int nNodes_owned,
                        bool simultaneous=True,
                        bool smoothBoundaries=True,
                        int[:] fixedNodesBoolArray=None,
                        double alpha=0.):
    """
    Centroid Smoothing:
    Mesh nodes are displaced to the centroid of the polygon/volume formed by
    neighbouring nodes

    Parameters
    ----------
    nodeArray_: double[:,:]
        array of mesh nodes
        (!) will be modified with smoothing
    nodeElementsOffsets: int[:]
        array to get offsets for neighbouring element numbers
    nodeElementsArray: int[:]
        array to get element number from nodeElementOffsets
    nodeMaterialTypes: int[:]
        array to know which node is interior node (materialType=0)
    elementBarycentersArray: double[:,:]
        barycenters of elements
    elementVolumesArray: double[:]
        volume of elements
    elementNodesArray: int[:,:]
        list of nodes per elements
    nNodes_owned: int
        number of nodes owned
    simultaneous: bool
        if True: simultaneous smoothing
        if False: sequential smoothing
                  (!) will update areas/volumes and barycenters in this case
    smoothBoundaries: bool
        if True: boundaries are smoothed (only with surrounding boundary nodes)
        if False: boudnary nodes are fixed
    fixedNodesBoolArray: int[:]
        array of same length as nodeArray with:
        fixedNodesBoolArray[node] = 1 -> node is fixed
        fixedNodesBoolArray[node] = 0 -> node is not fixed
    alpha: double
        penalty term related to original position of nodes (0 <= alpha <= 1)
    """
    cySmoothNodesCentroid(nodeArray_=nodeArray_,
                          nodeElementOffsets=nodeElementOffsets,
                          nodeElementsArray=nodeElementsArray,
                          nodeMaterialTypes=nodeMaterialTypes,
                          elementBarycentersArray=elementBarycentersArray,
                          elementVolumesArray=elementVolumesArray,
                          elementNodesArray=elementNodesArray,
                          nNodes_owned=nNodes_owned,
                          fixedNodesBoolArray=fixedNodesBoolArray,
                          simultaneous=simultaneous,
                          smoothBoundaries=smoothBoundaries,
                          alpha=alpha)

def updateDilationElements(double[:] elementDilationArray_,
                           double[:] elementVolumeArray,
                           double[:] elementVolumeTargetArray):
    cyUpdateDilationElements(elementDilationArray_=elementDilationArray_,
                             elementVolumeArray=elementVolumeArray,
                             elementVolumeTargetArray=elementVolumeTargetArray)

def getDilationElements(double[:] elementVolumeArray,
                        double[:] elementVolumeTargetArray):
    elementDilationArray_ = np.zeros(len(elementVolumeArray))
    cyUpdateDilationElements(elementDilationArray_=elementDilationArray_,
                             elementVolumeArray=elementVolumeArray,
                             elementVolumeTargetArray=elementVolumeTargetArray)
    return elementDilationArray_

def updateDistortionElements(double[:] elementDistortionArray_,
                             double[:,:,:,:] J_array,
                             double[:,:] detJ_array,
                             int nd):
    cyUpdateDistortionElements(elementDistortionArray_=elementDistortionArray_,
                               J_array=J_array,
                               detJ_array=detJ_array,
                               nd=nd)

def getDistortionElements(double[:,:,:,:] J_array,
                          double[:,:] detJ_array,
                          int nd):
    elementDistortionArray_ = np.zeros(len(detJ_array))
    cyUpdateDistortionElements(elementDistortionArray_=elementDistortionArray_,
                               J_array=J_array,
                               detJ_array=detJ_array,
                               nd=nd)
    return elementDistortionArray_

def updateInverseMeanRatioTriangleElements(double[:] IMRElementsArray_,
                                           double[:,:] nodeArray,
                                           int[:,:] elementNodesArray):
    cyUpdateInverseMeanRatioTriangleElements(IMRElementsArray_=IMRElementsArray_,
                                             nodeArray=nodeArray,
                                             elementNodesArray=elementNodesArray)

def getInverseMeanRatioTriangleElements(double[:,:] nodeArray,
                                        int[:,:] elementNodesArray):
    IMRElementsArray_ = np.zeros(len(elementNodesArray))
    cyUpdateInverseMeanRatioTriangleElements(IMRElementsArray_=IMRElementsArray_,
                                             nodeArray=nodeArray,
                                             elementNodesArray=elementNodesArray)
    return IMRElementsArray_

def updateInverseMeanRatioTriangleNodes(double[:] IMRNodesArray_,
                                        double[:,:] nodeArray,
                                        int[:,:] elementNodesArray,
                                        int[:] nodeElementOffsets,
                                        int[:] nodeElementsArray,
                                        bool el_average=False):
    cyUpdateInverseMeanRatioTriangleNodes(IMRNodesArray_=IMRNodesArray_,
                                          nodeArray=nodeArray,
                                          elementNodesArray=elementNodesArray,
                                          nodeElementOffsets=nodeElementOffsets,
                                          nodeElementsArray=nodeElementsArray,
                                          el_average=el_average)

def getInverseMeanRatioTriangleNodes(double[:,:] nodeArray,
                                     int[:,:] elementNodesArray,
                                     int[:] nodeElementOffsets,
                                     int[:] nodeElementsArray,
                                     bool el_average=False):
    IMRNodesArray_ = np.zeros(len(nodeArray))
    cyUpdateInverseMeanRatioTriangleNodes(IMRNodesArray_=IMRNodesArray_,
                                          nodeArray=nodeArray,
                                          elementNodesArray=elementNodesArray,
                                          nodeElementOffsets=nodeElementOffsets,
                                          nodeElementsArray=nodeElementsArray,
                                          el_average=el_average)
    return IMRNodesArray_

def getInverseMeanRatioSingleTriangle(int node0,
                                      double[:,:] nodeArray,
                                      int[:,:] elementNodesArray,
                                      int[:] nodeElementOffsets,
                                      int[:] nodeElementsArray):
    return cyGetInverseMeanRatioSingleTriangle(node0=node0,
                                               nodeArray=nodeArray,
                                               elementNodesArray=elementNodesArray,
                                               nodeElementOffsets=nodeElementOffsets,
                                               nodeElementsArray=nodeElementsArray)

def smoothNodesQuality(double[:] distortion,
                       double[:] dilation,
                       double[:,:] nodeArray,
                       int nNodes_owned,
                       int[:] nodeMaterialTypes,
                       int[:] nodeElementOffsets,
                       int[:] nodeElementsArray,
                       int[:] elementNodesArray,
                       bool apply_directly=False):
    assert 1>2, 'smoothNodesQuality is work in progress, do not use'
    return cySmoothNodesQuality(distortion=distortion,
                                dilation=dilation,
                                nodeArray=nodeArray,
                                nNodes_owned=nNodes_owned,
                                nodeMaterialTypes=nodeMaterialTypes,
                                nodeElementOffsets=nodeElementOffsets,
                                nodeElementsArray=nodeElementsArray,
                                elementNodesArray=elementNodesArray,
                                apply_directly=False)

def getLocalNearestNode(double[:] coords,
                        double[:,:] nodeArray,
                        int[:] nodeStarOffsets,
                        int[:] nodeStarArray,
                        int node):
    """Finds nearest node to coordinates (local)

    Parameters
    ----------
    coords: array_like
        coordinates from which to find nearest node
    nodeArray: array_like
        array of fluid mesh node coordinates
    nodeStarOffsets: array_like
        array of offsets from nodes (range)
    nodeStarArray: array_like
        array of neighbouring nodes
    node: int
        first guess for nearest node

    Returns
    -------
    node: int
        nearest node index
    """
    return pyxGetLocalNearestNode(coords,
                                  nodeArray,
                                  nodeStarOffsets,
                                  nodeStarArray,
                                  node)

def getLocalNearestElement(double[:] coords,
                           double[:,:] elementBarycentersArray,
                           int[:,:] elementNeighborsArray,
                           int eN):
    """Finds nearest element to coordinates (local)

    Parameters
    ----------
    coords: double[:]
        coordinates from which to find nearest element
    elementBarycentersArray: double[:,:]
        array of mesh cell barycenter coordinates
    elementNeighborsArray: double[:,:]
        array of element neighbors
    eN: int
        first guess for nearest element

    Returns
    -------
    eN: int
        nearest element index
    """
    return pyxGetLocalNearestElement(coords,
                                     elementBarycentersArray,
                                     elementNeighborsArray,
                                     eN)

def getLocalNearestElementAroundNode(double[:] coords,
                                     int[:] nodeElementOffsets,
                                     int[:] nodeElementsArray,
                                     double[:,:] elementBarycentersArray,
                                     int node):
    """Finds nearest neighbouring element of node to coordinates (local)

    Parameters
    ----------
    coords: double[:]
        coordinates from which to find nearest element
    nodeElementOffsets: int[:]
        element offsets from nodes
    nodeElementsArray: int[:]
        elements array from nodeElementOffsets
    elementBarycentersArray: int[:]
        array of mesh cell barycenter coordinates
    node: int
        node from which to search

    Returns
    -------
    eN: int
        nearest element index
    """
    return pyxGetLocalNearestElementAroundNode(coords=coords,
                                               nodeElementOffsets=nodeElementOffsets,
                                               nodeElementsArray=nodeElementsArray,
                                               elementBarycentersArray=elementBarycentersArray,
                                               node=node)

def getLocalNearestElementIntersection(double[:] coords,
                                       double[:,:,:] elementBoundaryNormalsArray,
                                       int[:,:] elementBoundariesArray,
                                       double[:,:] elementBoundaryBarycentersArray,
                                       int[:,:] elementNeighborsArray,
                                       double[:,:] elementBarycentersArray,
                                       int[:] exteriorElementBoundariesBoolArray,
                                       int eN):
    return pyxGetLocalNearestElementIntersection(coords,
                                           elementBoundaryNormalsArray,
                                           elementBoundariesArray,
                                           elementBoundaryBarycentersArray,
                                           elementNeighborsArray,
                                           elementBarycentersArray,
                                           exteriorElementBoundariesBoolArray,
                                           eN)

def getLocalElement(femSpace,
                    coords,
                    node):
    """Given coordinates and its nearest node, determine if it is on a
    local element.
    (!) old implementation -> slow

    Parameters
    ----------
    femSpace: object
        finite element space
    coords: array_like
        coordinates from which to element
    node: int
        nearest node index

    Returns
    -------
    eN: int or None
        local index of element (None if not found)
    """
    patchBoundaryNodes=set()
    checkedElements=[]
    # nodeElementOffsets give the indices to get the elements sharing the node
    #log Profiling.logEvent("Getting Local Element")
    if node+1 < len(femSpace.mesh.nodeElementOffsets):
        for eOffset in range(femSpace.mesh.nodeElementOffsets[node], femSpace.mesh.nodeElementOffsets[node + 1]):
            eN = femSpace.mesh.nodeElementsArray[eOffset]
            checkedElements.append(eN)
            # union of set
            patchBoundaryNodes|=set(femSpace.mesh.elementNodesArray[eN])
            # evaluate the inverse map for element eN (global to local)
            xi = femSpace.elementMaps.getInverseValue(eN, coords)
            #J = femSpace.elementMaps.getJacobianValues(eN, )
            # query whether xi lies within the reference element
            if femSpace.elementMaps.referenceElement.onElement(xi):
                return eN
    else:
        for eOffset in range(femSpace.mesh.nodeElementOffsets[node]):
            eN = femSpace.mesh.nodeElementsArray[eOffset]
            checkedElements.append(eN)
            # union of set
            patchBoundaryNodes|=set(femSpace.mesh.elementNodesArray[eN])
            # evaluate the inverse map for element eN (global to local)
            xi = femSpace.elementMaps.getInverseValue(eN, coords)
            #J = femSpace.elementMaps.getJacobianValues(eN, )
            # query whether xi lies within the reference element
            if femSpace.elementMaps.referenceElement.onElement(xi):
                return eN
    # extra loop if case coords is in neighbour element
    for node in patchBoundaryNodes:
        for eOffset in range(femSpace.mesh.nodeElementOffsets[node], femSpace.mesh.nodeElementOffsets[node + 1]):
            eN = femSpace.mesh.nodeElementsArray[eOffset]
            if eN not in checkedElements:
                checkedElements.append(eN)
                # evaluate the inverse map for element eN
                xi = femSpace.elementMaps.getInverseValue(eN, coords)
                # query whether xi lies within the reference element
                if femSpace.elementMaps.referenceElement.onElement(xi):
                    return eN
    # no elements found
    return None


def updateElementBoundaryNormalsTriangle(double[:,:,:] elementBoundaryNormalsArray_,
                                         double[:,:] nodeArray,
                                         int[:,:] elementBoundariesArray,
                                         int[:,:] elementBoundaryNodesArray,
                                         double[:,:] elementBoundaryBarycentersArray,
                                         double[:,:] elementBarycentersArray):
    pyxUpdateElementBoundaryNormalsTriangle(elementBoundaryNormalsArray_=elementBoundaryNormalsArray_,
                                              nodeArray=nodeArray,
                                              elementBoundariesArray=elementBoundariesArray,
                                              elementBoundaryNodesArray=elementBoundaryNodesArray,
                                              elementBoundaryBarycentersArray=elementBoundaryBarycentersArray,
                                              elementBarycentersArray=elementBarycentersArray)

def getElementBoundaryNormalsTriangle(double[:,:] nodeArray,
                                      int[:,:] elementBoundariesArray,
                                      int[:,:] elementBoundaryNodesArray,
                                      double[:,:] elementBoundaryBarycentersArray,
                                      double[:,:] elementBarycentersArray):
    elementBoundaryNormalsArray_ = np.zeros(elementBoundariesArray.shape[0], 3, 3)
    pyxUpdateElementBoundaryNormalsTriangle(elementBoundaryNormalsArray_=elementBoundaryNormalsArray_,
                                              nodeArray=nodeArray,
                                              elementBoundariesArray=elementBoundariesArray,
                                              elementBoundaryNodesArray=elementBoundaryNodesArray,
                                              elementBoundaryBarycentersArray=elementBoundaryBarycentersArray,
                                              elementBarycentersArray=elementBarycentersArray)
    return elementBoundaryNormalsArray_

def updateElementBoundaryNormalsTetra(double[:,:,:] elementBoundaryNormalsArray_,
                                      double[:,:] nodeArray,
                                      int[:,:] elementBoundariesArray,
                                      int[:,:] elementBoundaryNodesArray,
                                      double[:,:] elementBoundaryBarycentersArray,
                                      double[:,:] elementBarycentersArray):
    pyxUpdateElementBoundaryNormalsTetra(elementBoundaryNormalsArray_=elementBoundaryNormalsArray_,
                                         nodeArray=nodeArray,
                                         elementBoundariesArray=elementBoundariesArray,
                                         elementBoundaryNodesArray=elementBoundaryNodesArray,
                                         elementBoundaryBarycentersArray=elementBoundaryBarycentersArray,
                                         elementBarycentersArray=elementBarycentersArray)

def getElementBoundaryNormalsTetra(double[:,:] nodeArray,
                                   int[:,:] elementBoundariesArray,
                                   int[:,:] elementBoundaryNodesArray,
                                   double[:,:] elementBoundaryBarycentersArray,
                                   double[:,:] elementBarycentersArray):
    elementBoundaryNormalsArray_ = np.zeros(elementBoundariesArray.shape[0], 4, 3)
    pyxUpdateElementBoundaryNormalsTetra(elementBoundaryNormalsArray_=elementBoundaryNormalsArray_,
                                         nodeArray=nodeArray,
                                         elementBoundariesArray=elementBoundariesArray,
                                         elementBoundaryNodesArray=elementBoundaryNodesArray,
                                         elementBoundaryBarycentersArray=elementBoundaryBarycentersArray,
                                         elementBarycentersArray=elementBarycentersArray)
    return elementBoundaryNormalsArray_

def updateElementVolumesTriangle(double[:] elementVolumesArray_,
                                 int[:,:] elementNodesArray,
                                 double[:,:] nodeArray):
    cyUpdateElementVolumesTriangle(elementVolumesArray_=elementVolumesArray_,
                                   elementNodesArray=elementNodesArray,
                                   nodeArray=nodeArray)

def getElementVolumeTriangle(double[:] nA,
                             double[:] nB,
                             double[:] nC):
    return cyGetElementVolumeTriangle(nA=nA,
                                      nB=nB,
                                      nC=nC)

def updateElementVolumesTetra(double[:] elementVolumesArray_,
                              int[:,:] elementNodesArray,
                              double[:,:] nodeArray):
    cyUpdateElementVolumesTetra(elementVolumesArray_=elementVolumesArray_,
                                elementNodesArray=elementNodesArray,
                                nodeArray=nodeArray)

def updateElementBarycenters(double[:,:] elementBarycentersArray_,
                             int[:,:] elementNodesArray,
                             double[:,:] nodeArray):
    cyUpdateElementBarycenters(elementBarycentersArray_=elementBarycentersArray_,
                               elementNodesArray=elementNodesArray,
                               nodeArray=nodeArray)

def getCornerNodesTriangle(double[:,:] nodeArray,
                           int[:] nodeStarArray,
                           int[:] nodeStarOffsets,
                           int[:] nodeMaterialTypes,
                           int nNodes_owned):
    return cyGetCornerNodesTriangle(nodeArray=nodeArray,
                                    nodeStarArray=nodeStarArray,
                                    nodeStarOffsets=nodeStarOffsets,
                                    nodeMaterialTypes=nodeMaterialTypes,
                                    nNodes_owned=nNodes_owned)

### Cython implementation of functions above

cdef void cySmoothNodesLaplace(double[:,:] nodeArray_,
                               int[:] nodeStarOffsets,
                               int[:] nodeStarArray,
                               int[:] nodeMaterialTypes,
                               int nNodes_owned,
                               bool simultaneous=False,
                               bool smoothBoundaries=True,
                               int[:] fixedNodesBoolArray=None,
                               double alpha=0.):
    cdef double[:,:] nodeArray0
    if simultaneous or alpha!= 0:
        nodeArray0 = np.zeros_like(nodeArray_)
        nodeArray0[:] = nodeArray_
    cdef double[:] sum_star = np.zeros(3)
    cdef double[:] vec = np.zeros(3)
    cdef double[:] vec2 = np.zeros(3)
    cdef double dot
    cdef double vec_dist
    cdef int nNodeInStar
    cdef int nNodes = 0
    cdef bool fixed_node = False
    for node in range(nNodes_owned):
        sum_star[:] = 0.
        nNodes = 0
        if nodeMaterialTypes[node] == 0:
            for nOffset in range(nodeStarOffsets[node],
                                 nodeStarOffsets[node+1]):
                if simultaneous is True:
                    sum_star[0] += nodeArray0[nodeStarArray[nOffset], 0]
                    sum_star[1] += nodeArray0[nodeStarArray[nOffset], 1]
                    sum_star[2] += nodeArray0[nodeStarArray[nOffset], 2]
                else:
                    sum_star[0] += nodeArray_[nodeStarArray[nOffset], 0]
                    sum_star[1] += nodeArray_[nodeStarArray[nOffset], 1]
                    sum_star[2] += nodeArray_[nodeStarArray[nOffset], 2]
                nNodes += 1
        elif smoothBoundaries is True:
            # tridelat todo: works only in 2D here
            # smooth on boundary only unless it is a corner node
            fixed_node = False
            if fixedNodesBoolArray is not None:
                if fixedNodesBoolArray[node] == 1:
                    sum_star[0] = nodeArray0[node, 0]
                    sum_star[1] = nodeArray0[node, 1]
                    sum_star[2] = nodeArray0[node, 2]
                    nNodes = 1
                    fixed_node = True
            if fixed_node is False:
                for nOffset in range(nodeStarOffsets[node],
                                    nodeStarOffsets[node+1]):
                    if nodeMaterialTypes[nodeStarArray[nOffset]] != 0:
                        if simultaneous is True:
                            sum_star[0] += nodeArray0[nodeStarArray[nOffset], 0]
                            sum_star[1] += nodeArray0[nodeStarArray[nOffset], 1]
                            sum_star[2] += nodeArray0[nodeStarArray[nOffset], 2]
                        else:
                            sum_star[0] += nodeArray_[nodeStarArray[nOffset], 0]
                            sum_star[1] += nodeArray_[nodeStarArray[nOffset], 1]
                            sum_star[2] += nodeArray_[nodeStarArray[nOffset], 2]
                        nNodes += 1
        if alpha != 0:
            nodeArray_[node, 0] = alpha*nodeArray0[node, 0]+(1-alpha)*sum_star[0]/nNodes
            nodeArray_[node, 1] = alpha*nodeArray0[node, 0]+(1-alpha)*sum_star[1]/nNodes
            nodeArray_[node, 2] = alpha*nodeArray0[node, 0]+(1-alpha)*sum_star[2]/nNodes
        else:
            nodeArray_[node, 0] = sum_star[0]/nNodes
            nodeArray_[node, 1] = sum_star[1]/nNodes
            nodeArray_[node, 2] = sum_star[2]/nNodes

cdef void cySmoothNodesCentroid(double[:,:] nodeArray_,
                                int[:] nodeElementOffsets,
                                int[:] nodeElementsArray,
                                int[:] nodeMaterialTypes,
                                double[:] elementVolumesArray,
                                double[:,:] elementBarycentersArray,
                                int[:,:] elementNodesArray,
                                int nNodes_owned,
                                int[:] fixedNodesBoolArray,
                                bool simultaneous=False,
                                bool smoothBoundaries = True,
                                double alpha=0.):
    cdef double[:,:] nodeArray0
    if simultaneous or alpha != 0:
        nodeArray0 = np.zeros_like(nodeArray_)
        nodeArray0[:] = nodeArray_
    cdef double[:] sum_star = np.zeros(3)
    cdef double[:] bN = np.zeros(3)
    cdef int nNodeInStar
    cdef double areas = 0.
    cdef int nNodes_star = 0
    cdef double[:] nodeOffset0
    cdef double[:] nodeOffset1
    cdef double var = 0.
    cdef int eN
    cdef double[:] centroid_cell
    cdef double[:] nA
    cdef double[:] nB
    cdef double[:] nC
    cdef double[:] vec = np.zeros(3)
    cdef double[:] vec2 = np.zeros(3)
    cdef double dot
    cdef double vec_dist
    for node in range(nNodes_owned):
        sum_star[:] = 0.
        areas = 0.
        if nodeMaterialTypes[node] == 0:
            for eOffset in range(nodeElementOffsets[node],
                                    nodeElementOffsets[node+1]):
                eN = nodeElementsArray[eOffset]
                sum_star[0] += elementBarycentersArray[eN,0]*elementVolumesArray[eN]
                sum_star[1] += elementBarycentersArray[eN,1]*elementVolumesArray[eN]
                sum_star[2] += elementBarycentersArray[eN,2]*elementVolumesArray[eN]
                areas += elementVolumesArray[eN]
            if alpha != 0:
                nodeArray_[node, 0] = alpha*nodeArray0[node, 0]+(1-alpha)*sum_star[0]/areas
                nodeArray_[node, 1] = alpha*nodeArray0[node, 0]+(1-alpha)*sum_star[1]/areas
                nodeArray_[node, 2] = alpha*nodeArray0[node, 0]+(1-alpha)*sum_star[2]/areas
            else:
                nodeArray_[node, 0] = sum_star[0]/areas
                nodeArray_[node, 1] = sum_star[1]/areas
                nodeArray_[node, 2] = sum_star[2]/areas
            if not simultaneous:
                for eOffset in range(nodeElementOffsets[node],
                                     nodeElementOffsets[node+1]):
                    eN = nodeElementsArray[eOffset]
                    nA = nodeArray_[elementNodesArray[eN, 0]]
                    nB = nodeArray_[elementNodesArray[eN, 1]]
                    nC = nodeArray_[elementNodesArray[eN, 2]]
                    elementVolumesArray[eN] = getElementVolumeTriangle(nA, nB, nC)
                    elementBarycentersArray[eN, 0] = (nA[0]+nB[0]+nC[0])/3.
                    elementBarycentersArray[eN, 1] = (nA[1]+nB[1]+nC[1])/3.
                    elementBarycentersArray[eN, 2] = (nA[2]+nB[2]+nC[2])/3.

# cdef tuple cyGetDistortionDilation(double[:,:,:,:] J_array,
#                                    double[:,:] detJ_array,
#                                    double[:] target_area_array,
#                                    int nd):
#     #W = np.array([[1., 0.5],
#     #              [0., 0.866]]) # [0., np.sqrt(3)/2.]
#     cdef double[:,:] J
#     cdef double[:,:] JT
#     cdef double detJ
#     cdef double trJTJ
#     cdef double[:,:] JTJ
#     cdef double dilation
#     cdef double distortion
#     cdef double[:] distortion_array = np.zeros(len(target_area_array))
#     cdef double[:] dilation_array = np.zeros(len(target_area_array))
#     for eN in range(len(target_area_array)):
#         detJ = detJ_array[eN, 0]  # index 0 as it is the same for all quadrature points
#         dilation = target_area_array[eN]/detJ
#         if dilation < 1.:
#             dilation = 1/dilation
#         dilation_array[eN] = dilation-1
#         J = J_array[eN][0]
#         JT = np.transpose(J)
#         JTJ = np.zeros_like(J)
#         JTJ = np.dot(J, JT)
#         trJTJ = 0.
#         for i in range(len(J)):
#             trJTJ += JTJ[i,i]
#         distortion = (1./nd*trJTJ)**(nd/2.)/detJ
#         distortion = 1-1./distortion
#         distortion_array[eN] = distortion
#         dilation_array[eN] = dilation
#     return distortion_array, dilation_array

cdef void cyUpdateDilationElements(double[:] elementDilationArray_,
                                   double[:] elementVolumeArray,
                                   double[:] elementVolumeTargetArray):
    cdef double dilation
    for eN in range(len(elementDilationArray_)):
        dilation = elementVolumeArray[eN]/elementVolumeTargetArray[eN]
        if dilation < 1.:
            dilation = 1/dilation
        elementDilationArray_[eN] = dilation

cdef void cyUpdateDistortionElements(double[:] elementDistortionArray_,
                                     double[:,:,:,:] J_array,
                                     double[:,:] detJ_array,
                                     int nd):
    cdef double[:,:] J
    cdef double[:,:] JT
    cdef double detJ
    cdef double trJTJ
    cdef double[:,:] JTJ
    for eN in range(len(elementDistortionArray_)):
        detJ = detJ_array[eN, 0]  # index 0 as it is the same for all quadrature points
        J = J_array[eN][0]
        JT = np.transpose(J)
        JTJ = np.zeros_like(J)
        JTJ = np.dot(J, JT)
        trJTJ = 0.
        for i in range(len(J)):
            trJTJ += JTJ[i,i]
        elementDistortionArray_[eN] = (1./nd*trJTJ)**(nd/2.)/detJ

cdef void cyUpdateInverseMeanRatioTriangleNodes(double[:] IMRNodesArray_,
                                                double[:,:] nodeArray,
                                                int[:,:] elementNodesArray,
                                                int[:] nodeElementOffsets,
                                                int[:] nodeElementsArray,
                                                bool el_average=False):
    cdef double[:,:] W = np.array([[1., 0.5],
                                   [0., np.sqrt(3)/2.]])
    cdef double[:,:] A = np.zeros((2,2))
    cdef double[:,:] AW = np.zeros((2,2))
    cdef double[:] IMR_elements
    if el_average:
        IMR_elements = np.zeros(len(elementNodesArray))
    cdef int[:] nElements = np.zeros(len(nodeArray), dtype=np.int32)
    cdef double[:] vec_a
    cdef double[:] vec_b
    cdef double[:] vec_c
    cdef double IMR
    cdef int nEl = 0
    for eN in range(len(elementNodesArray)):
        for iN, node in enumerate(elementNodesArray[eN]):
            vec_a = nodeArray[elementNodesArray[eN, iN]]
            vec_b = nodeArray[elementNodesArray[eN, iN-1]]
            vec_c = nodeArray[elementNodesArray[eN, iN-2]]
            A[0,0] = vec_b[0]-vec_a[0]
            A[1,0] = vec_b[1]-vec_a[1]
            A[0,1] = vec_c[0]-vec_a[0]
            A[1,1] = vec_c[1]-vec_a[1]
            AW = np.dot(A, np.linalg.inv(W))
            IMR = (AW[0,0]**2+AW[0,1]**2+AW[1,0]**2+AW[1,1]**2)/(2*np.abs(AW[0,0]*AW[1,1]-AW[0,1]*AW[1,0]))
            IMRNodesArray_[node] += IMR
            nElements[node] += 1
            if el_average:
                IMR_elements[eN] += IMR
        if el_average:
            IMR_elements[eN] = IMR_elements[eN]/3.
    for node in range(len(IMRNodesArray_)):
        if not el_average:
            IMRNodesArray_[node] = IMRNodesArray_[node]/nElements[node]
        else:
            nEl = 0
            IMRNodesArray_[node] = 0.
            for eOffset in range(nodeElementOffsets[node],
                                 nodeElementOffsets[node+1]):
                eN = nodeElementsArray[eOffset]
                nEl += 1
                IMRNodesArray_[node] += IMR_elements[eN]
            IMRNodesArray_[node] = IMRNodesArray_[node]/nEl

cdef void cyUpdateInverseMeanRatioTriangleElements(double[:] IMRElementsArray_,
                                                   double[:,:] nodeArray,
                                                   int[:,:] elementNodesArray):
    cdef double[:,:] W = np.array([[1., 0.5],
                                   [0., np.sqrt(3)/2.]])
    cdef double[:,:] A = np.zeros((2,2))
    cdef double[:,:] AW = np.zeros((2,2))
    cdef double[:] vec_a
    cdef double[:] vec_b
    cdef double[:] vec_c
    cdef double IMR
    cdef int nEl = 0
    for eN in range(len(elementNodesArray)):
        for iN, node in enumerate(elementNodesArray[eN]):
            vec_a = nodeArray[elementNodesArray[eN, iN]]
            vec_b = nodeArray[elementNodesArray[eN, iN-1]]
            vec_c = nodeArray[elementNodesArray[eN, iN-2]]
            A[0,0] = vec_b[0]-vec_a[0]
            A[1,0] = vec_b[1]-vec_a[1]
            A[0,1] = vec_c[0]-vec_a[0]
            A[1,1] = vec_c[1]-vec_a[1]
            AW = np.dot(A, np.linalg.inv(W))
            IMR = (AW[0,0]**2+AW[0,1]**2+AW[1,0]**2+AW[1,1]**2)/(2*np.abs(AW[0,0]*AW[1,1]-AW[0,1]*AW[1,0]))
            IMRElementsArray_[eN] += IMR
        IMRElementsArray_[eN] = IMRElementsArray_[eN]/3.

cdef double cyGetInverseMeanRatioSingleTriangle(int node0,
                                                double[:,:] nodeArray,
                                                int[:,:] elementNodesArray,
                                                int[:] nodeElementOffsets,
                                                int[:] nodeElementsArray,
                                                bool el_average=False):
    cdef double[:,:] W = np.array([[1., 0.5],
                                   [0., np.sqrt(3)/2.]])
    cdef double[:,:] A = np.zeros((2,2))
    cdef double[:,:] AW = np.zeros((2,2))
    cdef double IMR_node = 0.
    cdef int nEl = 0
    cdef double[:] vec_a
    cdef double[:] vec_b
    cdef double[:] vec_c
    cdef double IMR
    for eOffset in range(nodeElementOffsets[node0],
                         nodeElementOffsets[node0+1]):
        eN = nodeElementsArray[eOffset]
        nEl += 1
        ################
        for iN, node in enumerate(elementNodesArray[eN]):
            if el_average:
                vec_a = nodeArray[elementNodesArray[eN, iN]]
                vec_b = nodeArray[elementNodesArray[eN, iN-1]]
                vec_c = nodeArray[elementNodesArray[eN, iN-2]]
                A[0,0] = vec_b[0]-vec_a[0]
                A[1,0] = vec_b[1]-vec_a[1]
                A[0,1] = vec_c[0]-vec_a[0]
                A[1,1] = vec_c[1]-vec_a[1]
                AW = np.dot(A, np.linalg.inv(W))
                IMR = (AW[0,0]**2+AW[0,1]**2+AW[1,0]**2+AW[1,1]**2)/(2*np.abs(AW[0,0]*AW[1,1]-AW[0,1]*AW[1,0]))
                IMR_node += IMR/3.  # /3. because 3 nodes in element
            else:
                if node == node0:
                    vec_a = nodeArray[elementNodesArray[eN, iN]]
                    vec_b = nodeArray[elementNodesArray[eN, iN-1]]
                    vec_c = nodeArray[elementNodesArray[eN, iN-2]]
                    A[0,0] = vec_b[0]-vec_a[0]
                    A[1,0] = vec_b[1]-vec_a[1]
                    A[0,1] = vec_c[0]-vec_a[0]
                    A[1,1] = vec_c[1]-vec_a[1]
                    AW = np.dot(A, np.linalg.inv(W))
                    IMR = (AW[0,0]**2+AW[0,1]**2+AW[1,0]**2+AW[1,1]**2)/(2*np.abs(AW[0,0]*AW[1,1]-AW[0,1]*AW[1,0]))
                    IMR_node += IMR
    IMR_node = IMR_node/nEl
    return IMR_node

cdef double[:,:] cySmoothNodesQuality(double[:] distortion,
                                      double[:] dilation,
                                      double[:,:] nodeArray,
                                      int nNodes_owned,
                                      int[:] nodeMaterialTypes,
                                      int[:] nodeElementOffsets,
                                      int[:] nodeElementsArray,
                                      int[:] elementNodesArray,
                                      bool apply_directly=False):
    cdef double[:,:] disp = np.zeros_like(nodeArray)
    disp[:] = nodeArray
    cdef double[:,:] nodeArrayMod
    if not apply_directly:
        nodeArrayMod = np.zeros_like(nodeArray)
        nodeArrayMod[:] = nodeArray
    else:
        nodeArrayMod = nodeArray
    cdef double[:] weighted_pos = np.zeros(3)
    cdef double weight = 0
    cdef double weights = 0
    for node in range(nNodes_owned):
        if nodeMaterialTypes[node] == 0:
            weights = 0
            weighted_pos[0] = 0
            weighted_pos[1] = 0
            weighted_pos[2] = 0
            for eOffset in range(nodeElementOffsets[node],
                                 nodeElementOffsets[node+1]):
                eN = nodeElementsArray[eOffset]
                for eNnode in elementNodesArray[eN]:
                    if eNnode != node:
                        weight = distortion[eN]
                        weighted_pos[0] += nodeArrayMod[eNnode, 0]*weight
                        weighted_pos[1] += nodeArrayMod[eNnode, 1]*weight
                        weighted_pos[2] += nodeArrayMod[eNnode, 2]*weight
                        weights += weight
                nodeArrayMod[node, 0] = weighted_pos[0]/weights
                nodeArrayMod[node, 1] = weighted_pos[1]/weights
                nodeArrayMod[node, 2] = weighted_pos[2]/weights
        disp[node, 0] = disp[node, 0]-nodeArray[node, 0]
        disp[node, 1] = disp[node, 1]-nodeArray[node, 1]
        disp[node, 2] = disp[node, 2]-nodeArray[node, 2]
    return disp

cdef double[:] recoveryAtNodes(double[:] variable,
                               double[:] nodeElementsArray,
                               double[:] nodeElementOffsets):
    """
    variable:
         Variable in element
    """
    cdef double[:] recovered_variable = np.zeros(len(nodeElementOffsets))
    cdef int nb_el
    cdef double var_av
    for node in range(len(nodeElementOffsets)):
        nb_el = 0
        grad_av = 0.
        for eOffset in range(nodeElementOffsets[node],
                             nodeElementOffsets[node+1]):
            nb_el += 1
            eN = nodeElementsArray[eOffset]
            var_av += variable[eN]
        recovered_variable[node] = var_av/nb_el
    return recovered_variable


# cdef gradientRecoveryAtNodes(double[:] nodeElementArray,
#                              double[:] nodeElementOffsets
#                              double[:,:,:] grad):
#     cdef double[:] recovered_grad = np.zeros(nodeElenentArray)
#     cdef double[:] grad_eN
#     cdef double grad_av
#     cdef int nb_el
#     for node in range(len(self.mesh.nodeArray)):
#         nb_el = 0
#         grad_av = 0.
#         for eOffset in range(self.mesh.nodeElementOffsets[node],
#                                 self.mesh.nodeElementOffsets[node+1]):
#             nb_el += 1
#             grad_eN = grad[nodeElementsArray[eOffset], 0]  # same value at all quad points
#             grad_av += grad_eN
#         grad_av /= nb_el
#         self.grads[node] = grad_av
#     self.model.grads = self.grads

cdef int pyxGetLocalNearestNode(double[:] coords,
                                double[:,:] nodeArray,
                                int[:] nodeStarOffsets,
                                int[:] nodeStarArray,
                                int node):
    # determine local nearest node distance
    cdef int nearest_node = node
    cdef int nearest_node0 = node
    cdef int nOffsets
    cdef double dist
    cdef double min_dist
    cdef double[:] node_coords
    cdef bool found_node = False
    node_coords = nodeArray[node]
    min_dist = (node_coords[0]-coords[0])*(node_coords[0]-coords[0])+\
               (node_coords[1]-coords[1])*(node_coords[1]-coords[1])+\
               (node_coords[2]-coords[2])*(node_coords[2]-coords[2])
    cdef int i = 0
    while found_node is False:
        nearest_node0 = nearest_node
        for nOffset in range(nodeStarOffsets[nearest_node0],
                             nodeStarOffsets[nearest_node0+1]):
            node = nodeStarArray[nOffset]
            node_coords = nodeArray[node]
            dist = (node_coords[0]-coords[0])*(node_coords[0]-coords[0])+\
                   (node_coords[1]-coords[1])*(node_coords[1]-coords[1])+\
                   (node_coords[2]-coords[2])*(node_coords[2]-coords[2])
            if dist < min_dist:
                min_dist = dist
                nearest_node = node
        if nearest_node0 == nearest_node:
            found_node = True
        i += 1
    return nearest_node


cdef int pyxGetLocalNearestElement(double[:] coords,
                                   double[:,:] elementBarycentersArray,
                                   int[:,:] elementNeighborsArray,
                                   int eN):
    # determine local nearest node distance
    cdef int nearest_eN = eN
    cdef int nearest_eN0 = eN
    cdef int nOffsets
    cdef double dist
    cdef double min_dist
    cdef double[:] eN_coords
    cdef bool found_eN = False
    eN_coords = elementBarycentersArray[eN]
    min_dist = (eN_coords[0]-coords[0])*(eN_coords[0]-coords[0])+\
               (eN_coords[1]-coords[1])*(eN_coords[1]-coords[1])+\
               (eN_coords[2]-coords[2])*(eN_coords[2]-coords[2])
    cdef int i = 0
    while found_eN is False:
        nearest_eN0 = nearest_eN
        for eN in elementNeighborsArray[nearest_eN0]:
            eN_coords = elementBarycentersArray[eN]
            dist = (eN_coords[0]-coords[0])*(eN_coords[0]-coords[0])+\
                   (eN_coords[1]-coords[1])*(eN_coords[1]-coords[1])+\
                   (eN_coords[2]-coords[2])*(eN_coords[2]-coords[2])
            if dist < min_dist:
                min_dist = dist
                nearest_eN = eN
        if nearest_eN0 == nearest_eN:
            found_eN = True
        i += 1
    return nearest_eN

cdef int pyxGetLocalNearestElementIntersection(double[:] coords,
                                               double[:,:,:] elementBoundaryNormalsArray,
                                               int[:,:] elementBoundariesArray,
                                               double[:,:] elementBoundaryBarycentersArray,
                                               int[:,:] elementNeighborsArray,
                                               double[:,:] elementBarycentersArray,
                                               int[:] exteriorElementBoundariesBoolArray,
                                               int eN,
                                               double tol=1e-10):
    # determine local nearest node distance
    cdef int nearest_eN = eN
    cdef int nearest_eN0 = eN
    cdef int nOffsets
    cdef double dist
    cdef double min_dist
    cdef double[:] eN_coords = np.zeros(3)
    cdef int maxit = 2*len(elementBoundaryNormalsArray)
    eN_coords[0] = elementBarycentersArray[eN, 0]
    eN_coords[1] = elementBarycentersArray[eN, 1]
    eN_coords[2] = elementBarycentersArray[eN, 2]
    min_dist = np.sqrt((eN_coords[0]-coords[0])*(eN_coords[0]-coords[0])+\
                       (eN_coords[1]-coords[1])*(eN_coords[1]-coords[1])+\
                       (eN_coords[2]-coords[2])*(eN_coords[2]-coords[2]))
    cdef double[:] direction = np.zeros(3)
    direction[0] = (coords[0]-eN_coords[0])/min_dist
    direction[1] = (coords[1]-eN_coords[1])/min_dist
    direction[2] = (coords[2]-eN_coords[2])/min_dist
    cdef bool found_eN = False
    cdef int i = 0
    cdef double[:] bound_bar  # barycenter of boundary
    cdef double alpha  # distance
    cdef int b_i  # boundary index element
    cdef int b_i_last = -999
    cdef double[:] normal = np.zeros(3)  # element boundary normal
    cdef double dot  # dot product result
    cdef double dot2  # dot product 2 result
    cdef int it = 0
    while found_eN is False and it < maxit:
        nearest_eN0 = nearest_eN
        alpha_min = 1e12
        for j, b_i in enumerate(elementBoundariesArray[nearest_eN0]):
            if b_i != b_i_last and exteriorElementBoundariesBoolArray[b_i] == 0.:
                normal[0] = elementBoundaryNormalsArray[nearest_eN0, j, 0]
                normal[1] = elementBoundaryNormalsArray[nearest_eN0, j, 1]
                normal[2] = elementBoundaryNormalsArray[nearest_eN0, j, 2]
                bound_bar = elementBoundaryBarycentersArray[b_i]
                dot = normal[0]*direction[0]+normal[1]*direction[1]+normal[2]*direction[2]
                if dot > 0.:
                    dot2 = (bound_bar[0]-eN_coords[0])*normal[0]+(bound_bar[1]-eN_coords[1])*normal[1]+(bound_bar[2]-eN_coords[2])*normal[2]
                    if dot2 >= 0.:
                        alpha = dot2/dot
                        if 0. < alpha < alpha_min:
                            alpha_min = alpha
                            nearest_eN = elementNeighborsArray[nearest_eN0, j]
                            b_i_last = b_i
        if nearest_eN != nearest_eN0:
            if min_dist-alpha_min > 0:
                eN_coords[0] += alpha_min*direction[0]
                eN_coords[1] += alpha_min*direction[1]
                eN_coords[2] += alpha_min*direction[2]
                min_dist -= alpha_min
            else:  # going too far
                nearest_eN = nearest_eN0
        if nearest_eN0 == nearest_eN:
            found_eN = True
        i += 1
        it += 1
    if it >= maxit:
        assert 1>2, 'could not find element! (element: '+str(eN)+')'
    return nearest_eN

cdef int pyxGetLocalNearestElementAroundNode(double[:] coords,
                                             int[:] nodeElementOffsets,
                                             int[:] nodeElementsArray,
                                             double[:,:] elementBarycentersArray,
                                             int node):
    cdef double dist
    cdef double min_dist
    cdef double[:] eN_coords
    cdef int i = 0
    for eN in nodeElementsArray[nodeElementOffsets[node]:nodeElementOffsets[node+1]]:
        eN_coords = elementBarycentersArray[eN]
        dist = (eN_coords[0]-coords[0])*(eN_coords[0]-coords[0])+\
               (eN_coords[1]-coords[1])*(eN_coords[1]-coords[1])+\
               (eN_coords[2]-coords[2])*(eN_coords[2]-coords[2])
        if i == 0 or dist < min_dist:
            min_dist = dist
            nearest_eN = eN
        i += 1
    # determine local nearest node distance
    return nearest_eN

cdef void pyxUpdateElementBoundaryNormalsTetra(double[:,:,:] elementBoundaryNormalsArray_,
                                               double[:,:] nodeArray,
                                               int[:,:] elementBoundariesArray,
                                               int[:,:] elementBoundaryNodesArray,
                                               double[:,:] elementBoundaryBarycentersArray,
                                               double[:,:] elementBarycentersArray):
    cdef double[:] normal_check = np.zeros(3)
    cdef double[:] U = np.zeros(3)
    cdef double[:] V = np.zeros(3)
    cdef double lengthn
    cdef int b_i
    cdef double[:] node0
    cdef double[:] node1
    cdef double[:] node2
    for i in range(len(elementBoundaryNormalsArray_)):
        for j in range(len(elementBoundaryNormalsArray_[i])):
            b_i = elementBoundariesArray[i, j]
            node0 = nodeArray[elementBoundaryNodesArray[elementBoundariesArray[i, j],0]]
            node1 = nodeArray[elementBoundaryNodesArray[elementBoundariesArray[i, j],1]]
            node2 = nodeArray[elementBoundaryNodesArray[elementBoundariesArray[i, j],2]]
            U[0] = node1[0]-node0[0]
            U[1] = node1[1]-node0[1]
            U[2] = node1[2]-node0[2]
            V[0] = node2[0]-node0[0]
            V[1] = node2[1]-node0[1]
            V[2] = node2[2]-node0[2]
            elementBoundaryNormalsArray_[i,j,0] = U[1]*V[2]-U[2]*V[1]
            elementBoundaryNormalsArray_[i,j,1] = U[2]*V[0]-U[0]*V[2]
            elementBoundaryNormalsArray_[i,j,2] = U[0]*V[1]-U[1]*V[0]
            lenghtn = np.sqrt(elementBoundaryNormalsArray_[i,j,0]**2+elementBoundaryNormalsArray_[i,j,1]**2+elementBoundaryNormalsArray_[i,j,2]**2)
            elementBoundaryNormalsArray_[i,j,0] /= lenghtn
            elementBoundaryNormalsArray_[i,j,1] /= lenghtn
            elementBoundaryNormalsArray_[i,j,2] /= lenghtn
            normal_check[0] = elementBoundaryBarycentersArray[b_i][0]-elementBarycentersArray[i][0]
            normal_check[1] = elementBoundaryBarycentersArray[b_i][1]-elementBarycentersArray[i][1]
            normal_check[2] = elementBoundaryBarycentersArray[b_i][2]-elementBarycentersArray[i][2]
            if np.dot(elementBoundaryNormalsArray_[i,j], normal_check) < 0:
                elementBoundaryNormalsArray_[i,j,0] = -elementBoundaryNormalsArray_[i,j,0]
                elementBoundaryNormalsArray_[i,j,1] = -elementBoundaryNormalsArray_[i,j,1]
                elementBoundaryNormalsArray_[i,j,2] = -elementBoundaryNormalsArray_[i,j,2]

cdef void pyxUpdateElementBoundaryNormalsTriangle(double[:,:,:] elementBoundaryNormalsArray_,
                                                  double[:,:] nodeArray,
                                                  int[:,:] elementBoundariesArray,
                                                  int[:,:] elementBoundaryNodesArray,
                                                  double[:,:] elementBoundaryBarycentersArray,
                                                  double[:,:] elementBarycentersArray):
    cdef double[:] normal_check = np.zeros(3)
    cdef double[:] U = np.zeros(3)
    cdef double lengthn
    cdef int b_i
    cdef double[:] node0
    cdef double[:] node1
    for i in range(len(elementBoundaryNormalsArray_)):
        for j in range(len(elementBoundaryNormalsArray_[i])):
            b_i = elementBoundariesArray[i, j]
            node0 = nodeArray[elementBoundaryNodesArray[elementBoundariesArray[i, j],0]]
            node1 = nodeArray[elementBoundaryNodesArray[elementBoundariesArray[i, j],1]]
            U[0] = node1[0]-node0[0]
            U[1] = node1[1]-node0[1]
            elementBoundaryNormalsArray_[i,j,0] = -U[1]
            elementBoundaryNormalsArray_[i,j,1] = U[0]
            lenghtn = np.sqrt(elementBoundaryNormalsArray_[i,j,0]**2+elementBoundaryNormalsArray_[i,j,1]**2)
            elementBoundaryNormalsArray_[i,j,0] /= lenghtn
            elementBoundaryNormalsArray_[i,j,1] /= lenghtn
            normal_check[0] = elementBoundaryBarycentersArray[b_i][0]-elementBarycentersArray[i][0]
            normal_check[1] = elementBoundaryBarycentersArray[b_i][1]-elementBarycentersArray[i][1]
            if np.dot(elementBoundaryNormalsArray_[i,j], normal_check) < 0:
                elementBoundaryNormalsArray_[i,j,0] = -elementBoundaryNormalsArray_[i,j,0]
                elementBoundaryNormalsArray_[i,j,1] = -elementBoundaryNormalsArray_[i,j,1]

cdef void cyUpdateElementVolumesTriangle(double[:] elementVolumesArray_,
                                         int[:,:] elementNodesArray,
                                         double[:,:] nodeArray):
    cdef double[:] nA
    cdef double[:] nB
    cdef double[:] nC
    cdef double base
    cdef double height
    cdef int nEl = len(elementVolumesArray_)
    cdef int eN
    for eN in range(nEl):
        nA = nodeArray[elementNodesArray[eN, 0]]
        nB = nodeArray[elementNodesArray[eN, 1]]
        nC = nodeArray[elementNodesArray[eN, 2]]
        base = np.sqrt((nB[1]-nA[1])**2+(nB[0]-nA[0])**2)
        height = np.abs((nB[1]-nA[1])*nC[0]-(nB[0]-nA[0])*nC[1]+nB[0]*nA[1]-nB[1]*nA[0])/base
        elementVolumesArray_[eN] = 0.5*base*height

cdef double cyGetElementVolumeTriangle(double[:] nA,
                                       double[:] nB,
                                       double[:] nC):
    cdef double base = np.sqrt((nB[1]-nA[1])**2+(nB[0]-nA[0])**2)
    cdef double height = np.abs((nB[1]-nA[1])*nC[0]-(nB[0]-nA[0])*nC[1]+nB[0]*nA[1]-nB[1]*nA[0])/base
    return 0.5*base*height

cdef void cyUpdateElementVolumesTetra(double[:] elementVolumesArray_,
                                      int[:,:] elementNodesArray,
                                      double[:,:] nodeArray):
    cdef double[:] nA
    cdef double[:] nB
    cdef double[:] nC
    cdef double[:] nD
    cdef double base_tri
    cdef double height_tri
    cdef double area_tri
    cdef double height_tetra
    cdef int nEl = len(elementVolumesArray_)
    cdef int eN
    for eN in range(nEl):
        nA = nodeArray[elementNodesArray[eN, 0]]
        nB = nodeArray[elementNodesArray[eN, 1]]
        nC = nodeArray[elementNodesArray[eN, 2]]
        nD = nodeArray[elementNodesArray[eN, 3]]
        base_tri = np.sqrt((nB[1]-nA[1])**2+(nB[0]-nA[0])**2)
        height_tri = np.abs((nB[1]-nA[1])*nC[0]-(nB[0]-nA[0])*nC[1]+nB[0]*nA[1]-nB[1]*nA[0])/base_tri
        area_tri = 0.5*base_tri*height_tri
        height_tetra = 0.
        elementVolumesArray_[eN] = 1./3.*area_tri*height_tetra

cdef void cyUpdateElementBarycenters(double[:,:] elementBarycentersArray_,
                                     int[:,:] elementNodesArray,
                                     double[:,:] nodeArray):
    cdef int nEl = len(elementBarycentersArray_)
    cdef int eN
    for eN in range(nEl):
        elementBarycentersArray_[eN, :] = 0.
        for iN, node in enumerate(elementNodesArray[eN]):
            elementBarycentersArray_[eN, 0] += nodeArray[node, 0]
            elementBarycentersArray_[eN, 1] += nodeArray[node, 1]
            elementBarycentersArray_[eN, 2] += nodeArray[node, 2]
        elementBarycentersArray_[eN, 0] = elementBarycentersArray_[eN, 0]/(iN+1)
        elementBarycentersArray_[eN, 1] = elementBarycentersArray_[eN, 1]/(iN+1)
        elementBarycentersArray_[eN, 2] = elementBarycentersArray_[eN, 2]/(iN+1)

cdef np.ndarray cyGetCornerNodesTriangle(double[:,:] nodeArray,
                                         int[:] nodeStarArray,
                                         int[:] nodeStarOffsets,
                                         int[:] nodeMaterialTypes,
                                         int nNodes_owned):
    cdef np.ndarray cornerNodesArray = np.array([], dtype=np.int32)
    cdef double[:] vec = np.zeros(3)
    cdef double[:] vec2 = np.zeros(3)
    cdef double vec_dist
    cdef double dot
    for node in range(nNodes_owned):
        if nodeMaterialTypes[node] != 0:
            vec[:] = 0.
            for nOffset in range(nodeStarOffsets[node],
                                 nodeStarOffsets[node+1]):
                if nodeMaterialTypes[nodeStarArray[nOffset]] != 0:
                    if vec[0] == 0. and vec[1] == 0. and vec[2] == 0.:
                        # initialize first vector
                        vec[0] = nodeArray[node, 0]-nodeArray[nodeStarArray[nOffset], 0]
                        vec[1] = nodeArray[node, 1]-nodeArray[nodeStarArray[nOffset], 1]
                        vec[2] = nodeArray[node, 2]-nodeArray[nodeStarArray[nOffset], 2]
                        vec_dist = np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2)
                        vec[0] = vec[0]/vec_dist
                        vec[1] = vec[1]/vec_dist
                        vec[2] = vec[2]/vec_dist
                    else:
                        vec2[0] = nodeArray[node, 0]-nodeArray[nodeStarArray[nOffset], 0]
                        vec2[1] = nodeArray[node, 1]-nodeArray[nodeStarArray[nOffset], 1]
                        vec2[2] = nodeArray[node, 2]-nodeArray[nodeStarArray[nOffset], 2]
                        vec_dist = np.sqrt(vec2[0]**2+vec2[1]**2+vec2[2]**2)
                        vec2[0] = vec2[0]/vec_dist
                        vec2[1] = vec2[1]/vec_dist
                        vec2[2] = vec2[2]/vec_dist
                        dot = vec[0]*vec2[0]+vec[1]*vec2[1]+vec[2]*vec2[2]
                        if dot == 1. or dot == -1.:
                            dot = 1
                        else:
                            cornerNodesArray = np.append(cornerNodesArray, node)
    return cornerNodesArray
