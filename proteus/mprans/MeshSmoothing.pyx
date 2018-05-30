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
                       bool apply_directly=False):
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
    return cySmoothNodesLaplace(nodeArray,
                                nodeStarOffsets,
                                nodeStarArray,
                                nodeMaterialTypes,
                                nNodes_owned,
                                nSmooth,
                                boundaryNormals,
                                fixedNodes,
                                apply_directly)


cdef double[:,:] cySmoothNodesLaplace(double[:,:] nodeArray,
                                      int[:] nodeStarOffsets,
                                      int[:] nodeStarArray,
                                      int[:] nodeMaterialTypes,
                                      int nNodes_owned,
                                      int nSmooth,
                                      double[:,:] boundaryNormals,
                                      int[:] fixedNodes,
                                      bool apply_directly=False):
    cdef double[:,:] disp = np.zeros_like(nodeArray)
    cdef double[:,:] nodeArrayMod
    if not apply_directly:
        nodeArrayMod = np.zeros_like(nodeArray)
        nodeArrayMod[:] = nodeArray
    else:
        nodeArrayMod = nodeArray
    cdef double[:] sum_star = np.zeros(3)
    cdef int nNodeInStar
    for i in range(nSmooth):
        for node in range(nNodes_owned):
            sum_star[:] = 0.
            if nodeMaterialTypes[node] == 0:
                for nOffset in range(nodeStarOffsets[node],
                                     nodeStarOffsets[node+1]):
                    sum_star[0] += nodeArrayMod[nodeStarArray[nOffset], 0]
                    sum_star[1] += nodeArrayMod[nodeStarArray[nOffset], 1]
                    sum_star[2] += nodeArrayMod[nodeStarArray[nOffset], 2]
                nNodeInStar = abs(nodeStarOffsets[node]-nodeStarOffsets[node+1])
                nodeArrayMod[node, 0] = sum_star[0]/nNodeInStar
                nodeArrayMod[node, 1] = sum_star[1]/nNodeInStar
                nodeArrayMod[node, 2] = sum_star[2]/nNodeInStar
            else:
                if boundaryNormals is not None:
                    fixed = False
                    if fixedNodes is not None:
                        if fixedNodes[nodeMaterialTypes[node]] == 1:
                            fixed = True
                    if fixed is False:
                        bN = boundaryNormals[nodeMaterialTypes[node]]
                        if not bN[0] == 0 and bN[1] == 0 and bN[2] == 0:
                            for nOffset in range(nodeStarOffsets[node],
                                                 nodeStarOffsets[node+1]):
                                sum_star[0] += nodeArrayMod[nodeStarArray[nOffset], 0]
                                sum_star[1] += nodeArrayMod[nodeStarArray[nOffset], 1]
                                sum_star[2] += nodeArrayMod[nodeStarArray[nOffset], 2]
                            nNodeInStar = abs(nodeStarOffsets[node]-nodeStarOffsets[node+1])
                            nodeArrayMod[node, 0] = sum_star[0]/nNodeInStar*(1-bN[0])
                            nodeArrayMod[node, 1] = sum_star[1]/nNodeInStar*(1-bN[1])
                            nodeArrayMod[node, 2] = sum_star[2]/nNodeInStar*(1-bN[2])
            if i == nSmooth-1:
                disp[node, 0] = nodeArrayMod[node, 0]-nodeArray[node, 0]
                disp[node, 1] = nodeArrayMod[node, 1]-nodeArray[node, 1]
                disp[node, 2] = nodeArrayMod[node, 2]-nodeArray[node, 2]
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

cdef tuple cyGetQualityMetrics(double[:,:,:] J_array,
                               double[:,:] detJ_array,
                               double[:] target_area_array,
                               int nd):
    #W = np.array([[1., 0.5],
    #              [0., 0.866]]) # [0., np.sqrt(3)/2.]
    cdef double[:] J
    cdef double[:] JT
    cdef double detJ
    cdef double trJTJ
    cdef double[:,:] JTJ
    cdef double dilation
    cdef double distortion
    cdef double[:] distortion_array
    cdef double[:] dilation_array
    for eN in range(len(target_area_array)):
        detJ = detJ_array[eN, 0]  # index 0 as it is the same for all quadrature points
        dilation = target_area_array[eN]/detJ  
        if dilation < 1.:
            dilation = 1/dilation
        dilation_array[eN] = dilation-1
        J = J_array[eN][0]
        JT = J.T
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


def getQualityMetrics(double[:,:,:] J_array,
                      double[:,:] detJ_array,
                      double[:] target_area_array,
                      int nd):
    return cyGetQualityMetrics(J_array=J_array,
                             detJ_array=detJ_array,
                             target_area_array=target_area_array,
                             nd=nd)


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
        disp[node, 0] = nodeArrayMod[node, 0]-nodeArray[node, 0]
        disp[node, 1] = nodeArrayMod[node, 1]-nodeArray[node, 1]
        disp[node, 2] = nodeArrayMod[node, 2]-nodeArray[node, 2]
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
