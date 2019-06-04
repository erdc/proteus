# A type of -*- python -*- file
from proteus import Comm as proteus_Comm
import numpy as np
cimport numpy as np
from mpi4py.MPI cimport (Comm,
                         Op)
from mpi4py.libmpi cimport MPI_Comm
cimport mesh
from proteus cimport cmeshTools
from proteus.partitioning cimport (c_partitionElements,
                                   c_partitionNodes,
                                   c_partitionNodesFromTetgenFiles,
                                   buildQuadraticSubdomain2GlobalMappings_1d,
                                   buildQuadraticSubdomain2GlobalMappings_2d,
                                   buildQuadraticSubdomain2GlobalMappings_3d,
                                   buildQuadraticCubeSubdomain2GlobalMappings_3d,
                                   buildDiscontinuousGalerkinSubdomain2GlobalMappings)

def partitionElements(Comm comm, int nLayersOfOverlap, cmeshTools.CMesh cmesh, cmeshTools.CMesh subdomain_cmesh):
    cmesh.mesh.subdomainp = &subdomain_cmesh.mesh
    c_partitionElements(comm.ob_mpi,
                        cmesh.mesh,
                        nLayersOfOverlap)
    return (
        np.asarray(<int[:(comm.size+1)]> cmesh.mesh.elementOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.mesh.subdomainp.nElements_global]> cmesh.mesh.elementNumbering_subdomain2global),
        np.asarray(<int[:(comm.size+1)]> cmesh.mesh.nodeOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.mesh.subdomainp.nNodes_global]> cmesh.mesh.nodeNumbering_subdomain2global),
        np.asarray(<int[:(comm.size+1)]> cmesh.mesh.elementBoundaryOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.mesh.subdomainp.nElementBoundaries_global]> cmesh.mesh.elementBoundaryNumbering_subdomain2global),
        np.asarray(<int[:(comm.size+1)]> cmesh.mesh.edgeOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.mesh.subdomainp.nEdges_global]> cmesh.mesh.edgeNumbering_subdomain2global)
    )

def partitionNodes(Comm comm, int nLayersOfOverlap, cmeshTools.CMesh cmesh, cmeshTools.CMesh subdomain_cmesh):
    cmesh.mesh.subdomainp = &subdomain_cmesh.mesh
    c_partitionNodes(comm.ob_mpi,
                     cmesh.mesh,
                     nLayersOfOverlap)
    return (
        np.asarray(<int[:comm.size+1]> cmesh.mesh.elementOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.mesh.subdomainp.nElements_global]> cmesh.mesh.elementNumbering_subdomain2global),
        np.asarray(<int[:comm.size+1]> cmesh.mesh.nodeOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.mesh.subdomainp.nNodes_global]> cmesh.mesh.nodeNumbering_subdomain2global),
        np.asarray(<int[:comm.size+1]> cmesh.mesh.elementBoundaryOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.mesh.subdomainp.nElementBoundaries_global]> cmesh.mesh.elementBoundaryNumbering_subdomain2global),
        np.asarray(<int[:comm.size+1]> cmesh.mesh.edgeOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.mesh.subdomainp.nEdges_global]> cmesh.mesh.edgeNumbering_subdomain2global)
    )

def convertPUMIPartitionToPython(Comm comm, cmeshTools.CMesh cmesh, cmeshTools.CMesh subdomain_cmesh):
#need to have Python handles for the cmesh arrays created from apf/PUMI mesh 
    cmesh.mesh.subdomainp = &subdomain_cmesh.mesh
    return (
        np.asarray(<int[:comm.size+1]> cmesh.mesh.elementOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.mesh.subdomainp.nElements_global]> cmesh.mesh.elementNumbering_subdomain2global),
        np.asarray(<int[:comm.size+1]> cmesh.mesh.nodeOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.mesh.subdomainp.nNodes_global]> cmesh.mesh.nodeNumbering_subdomain2global),
        np.asarray(<int[:comm.size+1]> cmesh.mesh.elementBoundaryOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.mesh.subdomainp.nElementBoundaries_global]> cmesh.mesh.elementBoundaryNumbering_subdomain2global),
        np.asarray(<int[:comm.size+1]> cmesh.mesh.edgeOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.mesh.subdomainp.nEdges_global]> cmesh.mesh.edgeNumbering_subdomain2global)
    )

def partitionNodesFromTetgenFiles(Comm comm, object filebase, int indexBase, int nLayersOfOverlap, cmeshTools.CMesh cmesh, cmeshTools.CMesh subdomain_cmesh):
    cmesh.mesh.subdomainp = &subdomain_cmesh.mesh
    if not isinstance(filebase, bytes):
        filebase = filebase.encode()
    c_partitionNodesFromTetgenFiles(comm.ob_mpi,
                                    <const char*>(<char*>filebase),
                                    indexBase,
                                    cmesh.mesh,
                                    nLayersOfOverlap)
    return (
        np.asarray(<int[:comm.size+1]> cmesh.mesh.elementOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.mesh.subdomainp.nElements_global]> cmesh.mesh.elementNumbering_subdomain2global),
        np.asarray(<int[:comm.size+1]> cmesh.mesh.nodeOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.mesh.subdomainp.nNodes_global]> cmesh.mesh.nodeNumbering_subdomain2global),
        np.asarray(<int[:comm.size+1]> cmesh.mesh.elementBoundaryOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.mesh.subdomainp.nElementBoundaries_global]> cmesh.mesh.elementBoundaryNumbering_subdomain2global),
        np.asarray(<int[:comm.size+1]> cmesh.mesh.edgeOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.mesh.subdomainp.nEdges_global]> cmesh.mesh.edgeNumbering_subdomain2global)
    )

def buildQuadraticLocal2GlobalMappings(Comm comm,
                                       int nSpace,
                                       cmeshTools.CMesh cmesh,
                                       cmeshTools.CMesh subdomain_cmesh,
                                       np.ndarray elementOffsets_subdomain_owned,
                                       np.ndarray nodeOffsets_subdomain_owned,
                                       np.ndarray elementBoundaryOffsets_subdomain_owned,
                                       np.ndarray edgeOffsets_subdomain_owned,
                                       np.ndarray elementNumbering_subdomain2global,
                                       np.ndarray nodeNumbering_subdomain2global,
                                       np.ndarray elementBoundaryNumbering_subdomain2global,
                                       np.ndarray edgeNumbering_subdomain2global,
                                       np.ndarray quadratic_dof_offsets_subdomain_owned,
                                       np.ndarray quadratic_subdomain_l2g,
                                       np.ndarray quadraticNumbering_subdomain2global,
                                       np.ndarray quadratic_lagrangeNodes):
    cdef int nDOF_all_processes=0
    cdef int nDOF_subdomain=0
    cdef int max_dof_neighbors=0
    if nSpace == 1:
        buildQuadraticSubdomain2GlobalMappings_1d(comm.ob_mpi,
                                                  cmesh.mesh,
                                                  <int*>(elementOffsets_subdomain_owned.data),
                                                  <int*>(nodeOffsets_subdomain_owned.data),
                                                  <int*>(elementNumbering_subdomain2global.data),
                                                  <int*>(nodeNumbering_subdomain2global.data),
                                                  nDOF_all_processes,
                                                  nDOF_subdomain,
                                                  max_dof_neighbors,
                                                  <int*>(quadratic_dof_offsets_subdomain_owned.data),
                                                  <int*>(quadratic_subdomain_l2g.data),
                                                  <int*>(quadraticNumbering_subdomain2global.data),
                                                  <double*>(quadratic_lagrangeNodes.data));
    elif nSpace == 2:
        buildQuadraticSubdomain2GlobalMappings_2d(comm.ob_mpi,
                                                  cmesh.mesh,
                                                  <int*>(elementBoundaryOffsets_subdomain_owned.data),
                                                  <int*>(nodeOffsets_subdomain_owned.data),
                                                  <int*>(elementBoundaryNumbering_subdomain2global.data),
                                                  <int*>(nodeNumbering_subdomain2global.data),
                                                  nDOF_all_processes,
                                                  nDOF_subdomain,
                                                  max_dof_neighbors,
                                                  <int*>(quadratic_dof_offsets_subdomain_owned.data),
                                                  <int*>(quadratic_subdomain_l2g.data),
                                                  <int*>(quadraticNumbering_subdomain2global.data),
                                                  <double*>(quadratic_lagrangeNodes.data))
    else:
        buildQuadraticSubdomain2GlobalMappings_3d(comm.ob_mpi,
                                                  cmesh.mesh,
                                                  <int*>(edgeOffsets_subdomain_owned.data),
                                                  <int*>(nodeOffsets_subdomain_owned.data),
                                                  <int*>(edgeNumbering_subdomain2global.data),
                                                  <int*>(nodeNumbering_subdomain2global.data),
                                                  nDOF_all_processes,
                                                  nDOF_subdomain,
                                                  max_dof_neighbors,
                                                  <int*>(quadratic_dof_offsets_subdomain_owned.data),
                                                  <int*>(quadratic_subdomain_l2g.data),
                                                  <int*>(quadraticNumbering_subdomain2global.data),
                                                  <double*>(quadratic_lagrangeNodes.data))
    return (nDOF_all_processes,
            nDOF_subdomain,
            max_dof_neighbors)
    
def buildQuadraticCubeLocal2GlobalMappings(Comm comm,
                                           int nSpace,
                                           cmeshTools.CMesh cmesh,
                                           cmeshTools.CMesh subdomain_cmesh,
                                           np.ndarray elementOffsets_subdomain_owned,
                                           np.ndarray nodeOffsets_subdomain_owned,
                                           np.ndarray elementBoundaryOffsets_subdomain_owned,
                                           np.ndarray edgeOffsets_subdomain_owned,
                                           np.ndarray elementNumbering_subdomain2global,
                                           np.ndarray nodeNumbering_subdomain2global,
                                           np.ndarray elementBoundaryNumbering_subdomain2global,
                                           np.ndarray edgeNumbering_subdomain2global,
                                           np.ndarray quadratic_dof_offsets_subdomain_owned,
                                           np.ndarray quadratic_subdomain_l2g,
                                           np.ndarray quadraticNumbering_subdomain2global,
                                           np.ndarray quadratic_lagrangeNodes):
    cdef int nDOF_all_processes=0
    cdef int nDOF_subdomain=0
    cdef int max_dof_neighbors=0
    if nSpace == 1:
        assert(False),"buildQuadraticCubeSubdomain2GlobalMappings_1d not implemented!!"
    elif nSpace == 2:
        assert(False),"buildQuadraticCubeSubdomain2GlobalMappings_2d not implemented!!"
    else:
        buildQuadraticCubeSubdomain2GlobalMappings_3d(comm.ob_mpi,
                                                      cmesh.mesh,
                                                      <int*>(edgeOffsets_subdomain_owned.data),
                                                      <int*>(nodeOffsets_subdomain_owned.data),
                                                      <int*>(edgeNumbering_subdomain2global.data),
                                                      <int*>(nodeNumbering_subdomain2global.data),
                                                      nDOF_all_processes,
                                                      nDOF_subdomain,
                                                      max_dof_neighbors,
                                                      <int*>(quadratic_dof_offsets_subdomain_owned.data),
                                                      <int*>(quadratic_subdomain_l2g.data),
                                                      <int*>(quadraticNumbering_subdomain2global.data),
                                                      <double*>(quadratic_lagrangeNodes.data))
    return (nDOF_all_processes,
            nDOF_subdomain,
            max_dof_neighbors)

def buildDiscontinuousGalerkinLocal2GlobalMappings(Comm comm,
                                                   int nDOF_element,
                                                   cmeshTools.CMesh cmesh,
                                                   cmeshTools.CMesh subdomain_cmesh,
                                                   np.ndarray elementOffsets_subdomain_owned,
                                                   np.ndarray elementNumbering_subdomain2global,
                                                   np.ndarray dg_dof_offsets_subdomain_owned,
                                                   np.ndarray dg_subdomain_l2g,
                                                   np.ndarray dgNumbering_subdomain2global):
    cdef int nDOF_all_processes=0
    cdef int nDOF_subdomain=0
    cdef int max_dof_neighbors=0
    buildDiscontinuousGalerkinSubdomain2GlobalMappings(comm.ob_mpi,
                                                       cmesh.mesh,
                                                       <int*>(elementOffsets_subdomain_owned.data),
                                                       <int*>(elementNumbering_subdomain2global.data),
                                                       nDOF_element,
                                                       nDOF_all_processes,
                                                       nDOF_subdomain,
                                                       max_dof_neighbors,
                                                       <int*>(dg_dof_offsets_subdomain_owned.data),
                                                       <int*>(dg_subdomain_l2g.data),
                                                       <int*>(dgNumbering_subdomain2global.data))
    return (nDOF_all_processes,
            nDOF_subdomain,
            max_dof_neighbors)
