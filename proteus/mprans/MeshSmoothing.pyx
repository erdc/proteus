cimport numpy as np
import numpy as np
from libcpp cimport bool

# cdef extern from "../CompKernel.h":
#     cdef cppclass CompKernelSpaceMapping

#   void calculateMapping_element(const int eN,
# 				                       const int k,
# 				                       double* mesh_dof,
# 				                       int* mesh_l2g,
# 				                       double* mesh_trial_ref,
# 				                       double* mesh_grad_trial_ref,
# 				                       double* jac,
# 				                       double& jacDet,
# 				                       double* jacInv,
# 				                       double& x,
# 				                       double& y,
# 				                       double& z);


def smoothNodesLaplace(double[:,:] nodeArray,
                       int[:] nodeStarOffsets,
                       int[:] nodeStarArray,
                       int[:] nodeMaterialTypes,
                       int nNodes_owned,
                       int nSmooth=1,
                       double[:,:] boundaryNormals=None,
                       int[:] fixedNodes=None,
                       bool apply_directly=False,
                       bool simultaneous=False):
    """
    Laplace Smoothing:
    Mesh nodes are displaced to the centroid of neighbouring nodes

    Parameters
    ----------
    nodeArray: double[:,:]
        array of mesh nodes
    nodeStarOffsets: int[:]
        array to get offsets for neighbouring node numbers
    nodeStarArray: int[:]
        array of neighbouring nodes
    nNodes_owned: int
        number of nodes owned
    nSmooth: int
        number of times Laplace smoothing should be executed
    """
    return cySmoothNodesLaplace(nodeArray=nodeArray,
                                nodeStarOffsets=nodeStarOffsets,
                                nodeStarArray=nodeStarArray,
                                nodeMaterialTypes=nodeMaterialTypes,
                                nNodes_owned=nNodes_owned,
                                nSmooth=nSmooth,
                                boundaryNormals=boundaryNormals,
                                fixedNodes=fixedNodes,
                                apply_directly=apply_directly,
                                simultaneous=simultaneous)


cdef double[:,:] cySmoothNodesLaplace(double[:,:] nodeArray,
                                      int[:] nodeStarOffsets,
                                      int[:] nodeStarArray,
                                      int[:] nodeMaterialTypes,
                                      int nNodes_owned,
                                      int nSmooth,
                                      double[:,:] boundaryNormals,
                                      int[:] fixedNodes,
                                      bool apply_directly=False,
                                      bool simultaneous=False):
    cdef double[:,:] disp = np.zeros_like(nodeArray)
    cdef double[:,:] nodeArray0 = np.zeros_like(nodeArray)
    nodeArray0[:] = nodeArray
    cdef double[:,:] nodeArraySimultaneous = np.zeros_like(nodeArray)
    nodeArraySimultaneous[:] = nodeArray
    cdef double[:,:] nodeArrayMod
    if not apply_directly:
        nodeArrayMod = np.zeros_like(nodeArray)
        nodeArrayMod[:] = nodeArray
    else:
        nodeArrayMod = nodeArray
    cdef double[:] sum_star = np.zeros(3)
    cdef double[:] bN = np.zeros(3)
    cdef int nNodeInStar
    for i in range(nSmooth):
        for node in range(nNodes_owned):
            sum_star[:] = 0.
            if nodeMaterialTypes[node] == 0:
                for nOffset in range(nodeStarOffsets[node],
                                     nodeStarOffsets[node+1]):
                    if simultaneous is True:
                        sum_star[0] += nodeArraySimultaneous[nodeStarArray[nOffset], 0]
                        sum_star[1] += nodeArraySimultaneous[nodeStarArray[nOffset], 1]
                        sum_star[2] += nodeArraySimultaneous[nodeStarArray[nOffset], 2]
                    else:
                        sum_star[0] += nodeArrayMod[nodeStarArray[nOffset], 0]
                        sum_star[1] += nodeArrayMod[nodeStarArray[nOffset], 1]
                        sum_star[2] += nodeArrayMod[nodeStarArray[nOffset], 2]
                nNodeInStar = abs(nodeStarOffsets[node]-nodeStarOffsets[node+1])
                nodeArrayMod[node, 0] = sum_star[0]/nNodeInStar
                nodeArrayMod[node, 1] = sum_star[1]/nNodeInStar
                nodeArrayMod[node, 2] = sum_star[2]/nNodeInStar
            # else:
            #     if boundaryNormals is not None:
            #         fixed = False
            #         if fixedNodes is not None:
            #             if fixedNodes[nodeMaterialTypes[node]] == 1:
            #                 fixed = True
            #         if fixed is False:
            #             bN = boundaryNormals[nodeMaterialTypes[node]]
            #             if not bN[0] == 0 and bN[1] == 0 and bN[2] == 0:
            #                 for nOffset in range(nodeStarOffsets[node],
            #                                      nodeStarOffsets[node+1]):
            #                     sum_star[0] += nodeArrayMod[nodeStarArray[nOffset], 0]
            #                     sum_star[1] += nodeArrayMod[nodeStarArray[nOffset], 1]
            #                     sum_star[2] += nodeArrayMod[nodeStarArray[nOffset], 2]
            #                 nNodeInStar = abs(nodeStarOffsets[node]-nodeStarOffsets[node+1])
            #                 nodeArrayMod[node, 0] = sum_star[0]/nNodeInStar*(1-np.abs(bN[0]))
            #                 nodeArrayMod[node, 1] = sum_star[1]/nNodeInStar*(1-np.abs(bN[1]))
            #                 nodeArrayMod[node, 2] = sum_star[2]/nNodeInStar*(1-np.abs(bN[2]))
            if i == nSmooth-1:
                disp[node, 0] = nodeArrayMod[node, 0]-nodeArray0[node, 0]
                disp[node, 1] = nodeArrayMod[node, 1]-nodeArray0[node, 1]
                disp[node, 2] = nodeArrayMod[node, 2]-nodeArray0[node, 2]
        if simultaneous is True:
            # update array for next smoothing iteration
            nodeArraySimultaneous[:] = nodeArrayMod
    return disp




# cdef cyCalcuteMapping_element(gg)
# ck.calculateMapping_element(eN,
#                             k,
#                             mesh_dof,
#                             mesh_l2g,
#                             mesh_trial_ref,
#                             mesh_grad_trial_ref,
#                             jac,
#                             jacDet,
#                             jacInv,
#                             x,y,z);

cdef tuple cyGetQualityMetrics(double[:,:,:,:] J_array,
                               double[:,:] detJ_array,
                               double[:] target_area_array,
                               int nd):
    #W = np.array([[1., 0.5],
    #              [0., 0.866]]) # [0., np.sqrt(3)/2.]
    cdef double[:,:] J
    cdef double[:,:] JT
    cdef double detJ
    cdef double trJTJ
    cdef double[:,:] JTJ
    cdef double dilation
    cdef double distortion
    cdef double[:] distortion_array = np.zeros(len(target_area_array))
    cdef double[:] dilation_array = np.zeros(len(target_area_array))
    for eN in range(len(target_area_array)):
        detJ = detJ_array[eN, 0]  # index 0 as it is the same for all quadrature points
        dilation = target_area_array[eN]/detJ
        if dilation < 1.:
            dilation = 1/dilation
        dilation_array[eN] = dilation-1
        J = J_array[eN][0]
        JT = np.transpose(J)
        JTJ = np.zeros_like(J)
        JTJ = np.dot(J, JT)
        trJTJ = 0.
        for i in range(len(J)):
            trJTJ += JTJ[i,i]
        distortion = (1./nd*trJTJ)**(nd/2.)/detJ
        distortion = 1-1./distortion
        distortion_array[eN] = distortion
        dilation_array[eN] = dilation
    return distortion_array, dilation_array

def getQualityMetrics(double[:,:,:,:] J_array,
                      double[:,:] detJ_array,
                      double[:] target_area_array,
                      int nd):
    return cyGetQualityMetrics(J_array=J_array,
                               detJ_array=detJ_array,
                               target_area_array=target_area_array,
                               nd=nd)

cdef tuple cyGetInverseMeanRatioTriangle(double[:,:] nodeArray,
                                         int[:,:] elementNodesArray,
                                         int[:] nodeElementOffsets,
                                         int[:] nodeElementsArray,
                                         bool el_average=False):
    cdef double[:,:] W = np.array([[1., 0.5],
                                   [0., np.sqrt(3)/2.]])
    cdef double[:,:] A = np.zeros((2,2))
    cdef double[:,:] AW = np.zeros((2,2))
    cdef double[:] IMR_nodes = np.zeros(len(nodeArray))
    cdef double[:] IMR_elements = np.zeros(len(elementNodesArray))
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
            IMR_nodes[node] += IMR
            nElements[node] += 1
            IMR_elements[eN] += IMR
        IMR_elements[eN] = IMR_elements[eN]/3.
    for node in range(len(IMR_nodes)):
        if not el_average:
            IMR_nodes[node] = IMR_nodes[node]/nElements[node]
        else:
            nEl = 0
            IMR_nodes[node] = 0.
            for eOffset in range(nodeElementOffsets[node],
                                nodeElementOffsets[node+1]):
                eN = nodeElementsArray[eOffset]
                nEl += 1
                IMR_nodes[node] += IMR_elements[eN]
            IMR_nodes[node] = IMR_nodes[node]/nEl
    return IMR_nodes, IMR_elements

def getInverseMeanRatioTriangle(double[:,:] nodeArray,
                                int[:,:] elementNodesArray,
                                int[:] nodeElementOffsets,
                                int[:] nodeElementsArray):
    return cyGetInverseMeanRatioTriangle(nodeArray=nodeArray,
                                         elementNodesArray=elementNodesArray,
                                         nodeElementOffsets=nodeElementOffsets,
                                         nodeElementsArray=nodeElementsArray)

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

def smoothNodesQuality(double[:] distortion,
                       double[:] dilation,
                       double[:,:] nodeArray,
                       int nNodes_owned,
                       int[:] nodeMaterialTypes,
                       int[:] nodeElementOffsets,
                       int[:] nodeElementsArray,
                       int[:] elementNodesArray,
                       bool apply_directly=False):
    return cySmoothNodesQuality(distortion=distortion,
                                dilation=dilation,
                                nodeArray=nodeArray,
                                nNodes_owned=nNodes_owned,
                                nodeMaterialTypes=nodeMaterialTypes,
                                nodeElementOffsets=nodeElementOffsets,
                                nodeElementsArray=nodeElementsArray,
                                elementNodesArray=elementNodesArray,
                                apply_directly=False)


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


cdef tuple pyxGetLocalNearestNode(double[:] coords,
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
    dist: float
        distance to nearest node
    """
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
    return nearest_node, dist

def getLocalNearestNode(double[:] coords,
                        double[:,:] nodeArray,
                        int[:] nodeStarOffsets,
                        int[:] nodeStarArray,
                        int node):
    return pyxGetLocalNearestNode(coords,
                                  nodeArray,
                                  nodeStarOffsets,
                                  nodeStarArray,
                                  node)


cdef tuple pyxGetLocalNearestElement(double[:] coords,
                                     double[:,:] elementBarycentersArray,
                                     int[:,:] elementNeighborsArray,
                                     int eN):
    """Finds nearest element to coordinates (local)
    Parameters
    ----------
    coords: array_like
        coordinates from which to find nearest element
    elementBarycentersArray: array_like
        array of mesh cell barycenter coordinates
    elementNeighborsArray: array_like
        array of element neighbors
    eN: int
        first guess for nearest element

    Returns
    -------
    eN: int
        nearest element index
    dist: float
        distance to nearest element
    """
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
    return nearest_eN, dist

def getLocalNearestElement(double[:] coords,
                           double[:,:] elementBarycentersArray,
                           int[:,:] elementNeighborsArray,
                           int eN):
    return pyxGetLocalNearestElement(coords,
                                     elementBarycentersArray,
                                     elementNeighborsArray,
                                     eN)

cdef tuple pyxGetLocalNearestElementIntersection(double[:] coords,
                                                 double[:,:,:] elementNormalsArray,
                                                 int[:,:] elementBoundariesArray,
                                                 double[:,:] elementBoundaryBarycentersArray,
                                                 int[:,:] elementNeighborsArray,
                                                 double[:,:] elementBarycentersArray,
                                                 int[:] exteriorElementBoundariesBoolArray,
                                                 int eN,
                                                 double tol=1e-10):
    """Finds nearest node to coordinates (local)
    """
    # determine local nearest node distance
    cdef int nearest_eN = eN
    cdef int nearest_eN0 = eN
    cdef int nOffsets
    cdef double dist
    cdef double min_dist
    cdef double[:] eN_coords = np.zeros(3)
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
    while found_eN is False:
        nearest_eN0 = nearest_eN
        alpha_min = 1e12
        for j, b_i in enumerate(elementBoundariesArray[nearest_eN0]):
            if b_i != b_i_last and exteriorElementBoundariesBoolArray[b_i] == 0.:
                normal[0] = elementNormalsArray[nearest_eN0, j, 0]
                normal[1] = elementNormalsArray[nearest_eN0, j, 1]
                normal[2] = elementNormalsArray[nearest_eN0, j, 2]
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
    return nearest_eN, min_dist

def getLocalNearestElementIntersection(double[:] coords,
                                       double[:,:,:] elementNormalsArray,
                                       int[:,:] elementBoundariesArray,
                                       double[:,:] elementBoundaryBarycentersArray,
                                       int[:,:] elementNeighborsArray,
                                       double[:,:] elementBarycentersArray,
                                       int[:] exteriorElementBoundariesBoolArray,
                                       int eN):
    return pyxGetLocalNearestElementIntersection(coords,
                                           elementNormalsArray,
                                           elementBoundariesArray,
                                           elementBoundaryBarycentersArray,
                                           elementNeighborsArray,
                                           elementBarycentersArray,
                                           exteriorElementBoundariesBoolArray,
                                           eN)

cdef tuple pyxGetLocalNearestElementAroundNode(double[:] coords,
                                               int[:] nodeElementOffsets,
                                               int[:] nodeElementsArray,
                                               double[:,:] elementBarycentersArray,
                                               int node):
    """
    """
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
    return nearest_eN, min_dist

def getLocalNearestElementAroundNode(coords,
                                     nodeElementOffsets,
                                     nodeElementsArray,
                                     elementBarycentersArray,
                                     node):
    return pyxGetLocalNearestElementAroundNode(coords=coords,
                                               nodeElementOffsets=nodeElementOffsets,
                                               nodeElementsArray=nodeElementsArray,
                                               elementBarycentersArray=elementBarycentersArray,
                                               node=node)

def getLocalElement(femSpace, coords, node):
    """Given coordinates and its nearest node, determine if it is on a
    local element.

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

cdef double[:,:,:] pyxGetElementBoundaryNormalsTetra3D(double[:,:] nodeArray,
                                                       int[:,:] elementBoundariesArray,
                                                       int[:,:] elementBoundaryNodesArray,
                                                       double[:,:] elementBoundaryBarycentersArray,
                                                       double[:,:] elementBarycentersArray):
    cdef double[:,:,:] normals = np.zeros((len(elementBoundariesArray), 4, 3))
    cdef double[:] normal_check = np.zeros(3)
    cdef double[:] U = np.zeros(3)
    cdef double[:] V = np.zeros(3)
    cdef double lengthn
    cdef int b_i
    cdef double[:] node0
    cdef double[:] node1
    cdef double[:] node2
    for i in range(len(normals)):
        for j in range(len(normals[i])):
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
            normals[i,j,0] = U[1]*V[2]-U[2]*V[1]
            normals[i,j,1] = U[2]*V[0]-U[0]*V[2]
            normals[i,j,2] = U[0]*V[1]-U[1]*V[0]
            lenghtn = np.sqrt(normals[i,j,0]**2+normals[i,j,1]**2+normals[i,j,2]**2)
            normals[i,j,0] /= lenghtn
            normals[i,j,1] /= lenghtn
            normals[i,j,2] /= lenghtn
            normal_check[0] = elementBoundaryBarycentersArray[b_i][0]-elementBarycentersArray[i][0]
            normal_check[1] = elementBoundaryBarycentersArray[b_i][1]-elementBarycentersArray[i][1]
            normal_check[2] = elementBoundaryBarycentersArray[b_i][2]-elementBarycentersArray[i][2]
            if np.dot(normals[i,j], normal_check) < 0:
                normals[i,j,0] = -normals[i,j,0]
                normals[i,j,1] = -normals[i,j,1]
                normals[i,j,2] = -normals[i,j,2]
    return normals

def getElementBoundaryNormalsTetra3D(nodeArray,
                                     elementBoundariesArray,
                                     elementBoundaryNodesArray,
                                     elementBoundaryBarycentersArray,
                                     elementBarycentersArray):
    return pyxGetElementBoundaryNormalsTetra3D(nodeArray=nodeArray,
                                               elementBoundariesArray=elementBoundariesArray,
                                               elementBoundaryNodesArray=elementBoundaryNodesArray,
                                               elementBoundaryBarycentersArray=elementBoundaryBarycentersArray,
                                               elementBarycentersArray=elementBarycentersArray)

cdef double[:,:,:] pyxGetElementBoundaryNormalsTriangle2D(double[:,:] nodeArray,
                                                          int[:,:] elementBoundariesArray,
                                                          int[:,:] elementBoundaryNodesArray,
                                                          double[:,:] elementBoundaryBarycentersArray,
                                                          double[:,:] elementBarycentersArray):
    cdef double[:,:,:] normals = np.zeros((len(elementBoundariesArray), 3, 3))
    cdef double[:] normal_check = np.zeros(3)
    cdef double[:] U = np.zeros(3)
    cdef double lengthn
    cdef int b_i
    cdef double[:] node0
    cdef double[:] node1
    for i in range(len(normals)):
        for j in range(len(normals[i])):
            b_i = elementBoundariesArray[i, j]
            node0 = nodeArray[elementBoundaryNodesArray[elementBoundariesArray[i, j],0]]
            node1 = nodeArray[elementBoundaryNodesArray[elementBoundariesArray[i, j],1]]
            U[0] = node1[0]-node0[0]
            U[1] = node1[1]-node0[1]
            normals[i,j,0] = -U[1]
            normals[i,j,1] = U[0]
            lenghtn = np.sqrt(normals[i,j,0]**2+normals[i,j,1]**2)
            normals[i,j,0] /= lenghtn
            normals[i,j,1] /= lenghtn
            normal_check[0] = elementBoundaryBarycentersArray[b_i][0]-elementBarycentersArray[i][0]
            normal_check[1] = elementBoundaryBarycentersArray[b_i][1]-elementBarycentersArray[i][1]
            if np.dot(normals[i,j], normal_check) < 0:
                normals[i,j,0] = -normals[i,j,0]
                normals[i,j,1] = -normals[i,j,1]
    return normals

def updateElementBoundaryNormalsTriangle2D(normals,
                                           nodeArray,
                                           elementBoundariesArray,
                                           elementBoundaryNodesArray,
                                           elementBoundaryBarycentersArray,
                                           elementBarycentersArray):
    return pyxUpdateElementBoundaryNormalsTriangle2D(normals=normals,
                                                     nodeArray=nodeArray,
                                                     elementBoundariesArray=elementBoundariesArray,
                                                     elementBoundaryNodesArray=elementBoundaryNodesArray,
                                                     elementBoundaryBarycentersArray=elementBoundaryBarycentersArray,
                                                     elementBarycentersArray=elementBarycentersArray)


cdef pyxUpdateElementBoundaryNormalsTetra3D(double[:,:,:] normals,
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
    for i in range(len(normals)):
        for j in range(len(normals[i])):
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
            normals[i,j,0] = U[1]*V[2]-U[2]*V[1]
            normals[i,j,1] = U[2]*V[0]-U[0]*V[2]
            normals[i,j,2] = U[0]*V[1]-U[1]*V[0]
            lenghtn = np.sqrt(normals[i,j,0]**2+normals[i,j,1]**2+normals[i,j,2]**2)
            normals[i,j,0] /= lenghtn
            normals[i,j,1] /= lenghtn
            normals[i,j,2] /= lenghtn
            normal_check[0] = elementBoundaryBarycentersArray[b_i][0]-elementBarycentersArray[i][0]
            normal_check[1] = elementBoundaryBarycentersArray[b_i][1]-elementBarycentersArray[i][1]
            normal_check[2] = elementBoundaryBarycentersArray[b_i][2]-elementBarycentersArray[i][2]
            if np.dot(normals[i,j], normal_check) < 0:
                normals[i,j,0] = -normals[i,j,0]
                normals[i,j,1] = -normals[i,j,1]
                normals[i,j,2] = -normals[i,j,2]

def updateElementBoundaryNormalsTetra3D(normals,
                                        nodeArray,
                                        elementBoundariesArray,
                                        elementBoundaryNodesArray,
                                        elementBoundaryBarycentersArray,
                                        elementBarycentersArray):
    return pyxUpdateElementBoundaryNormalsTetra3D(normals=normals,
                                                  nodeArray=nodeArray,
                                                  elementBoundariesArray=elementBoundariesArray,
                                                  elementBoundaryNodesArray=elementBoundaryNodesArray,
                                                  elementBoundaryBarycentersArray=elementBoundaryBarycentersArray,
                                                  elementBarycentersArray=elementBarycentersArray)

cdef pyxUpdateElementBoundaryNormalsTriangle2D(double[:,:,:] normals,
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
    for i in range(len(normals)):
        for j in range(len(normals[i])):
            b_i = elementBoundariesArray[i, j]
            node0 = nodeArray[elementBoundaryNodesArray[elementBoundariesArray[i, j],0]]
            node1 = nodeArray[elementBoundaryNodesArray[elementBoundariesArray[i, j],1]]
            U[0] = node1[0]-node0[0]
            U[1] = node1[1]-node0[1]
            normals[i,j,0] = -U[1]
            normals[i,j,1] = U[0]
            lenghtn = np.sqrt(normals[i,j,0]**2+normals[i,j,1]**2)
            normals[i,j,0] /= lenghtn
            normals[i,j,1] /= lenghtn
            normal_check[0] = elementBoundaryBarycentersArray[b_i][0]-elementBarycentersArray[i][0]
            normal_check[1] = elementBoundaryBarycentersArray[b_i][1]-elementBarycentersArray[i][1]
            if np.dot(normals[i,j], normal_check) < 0:
                normals[i,j,0] = -normals[i,j,0]
                normals[i,j,1] = -normals[i,j,1]

def updateElementBoundaryNormalsTriangle2D(normals,
                                           nodeArray,
                                           elementBoundariesArray,
                                           elementBoundaryNodesArray,
                                           elementBoundaryBarycentersArray,
                                           elementBarycentersArray):
    return pyxUpdateElementBoundaryNormalsTriangle2D(normals=normals,
                                                     nodeArray=nodeArray,
                                                     elementBoundariesArray=elementBoundariesArray,
                                                     elementBoundaryNodesArray=elementBoundaryNodesArray,
                                                     elementBoundaryBarycentersArray=elementBoundaryBarycentersArray,
                                                     elementBarycentersArray=elementBarycentersArray)
