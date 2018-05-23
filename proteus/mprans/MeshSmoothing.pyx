cimport numpy as np
import numpy as np

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


def smoothLaplace(double[:,:] nodeArray,
                  int[:] nodeStarOffsets,
                  int[:] nodeStarArray,
                  int[:] nodeMaterialTypes,
                  int nNodes_owned,
                  int nSmooth=1):
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
    cySmoothLaplace(nodeArray,
                    nodeStarOffsets,
                    nodeStarArray,
                    nodeMaterialTypes,
                    nNodes_owned,
                    nSmooth)


cdef cySmoothLaplace(double[:,:] nodeArray,
                     int[:] nodeStarOffsets,
                     int[:] nodeStarArray,
                     int[:] nodeMaterialTypes,
                     int nNodes_owned,
                     int nSmooth):
    cdef double[:] sum_star = np.zeros(3)
    cdef int nNodeInStar
    for i in range(nSmooth):
        for node in range(nNodes_owned):
            sum_star[:] = 0.
            if nodeMaterialTypes[node] == 0:
                for nOffset in range(nodeStarOffsets[node],
                                     nodeStarOffsets[node+1]):
                    sum_star[0] += nodeArray[nodeStarArray[nOffset]][0]
                    sum_star[1] += nodeArray[nodeStarArray[nOffset]][1]
                    sum_star[2] += nodeArray[nodeStarArray[nOffset]][2]
                nNodeInStar = abs(nodeStarOffsets[node]-nodeStarOffsets[node+1])
                nodeArray[node][0] = sum_star[0]/nNodeInStar
                nodeArray[node][1] = sum_star[1]/nNodeInStar
                nodeArray[node][2] = sum_star[2]/nNodeInStar

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
