# A type of -*- python -*- file
import Comm
import numpy as np
cimport numpy as np
from mpi4py.MPI cimport Comm
from mpi4py.libmpi cimport MPI_Comm
cimport mesh
from proteus cimport cmeshToolsTest

from proteus.partitioning cimport (c_partitionElements,
                                   c_partitionNodes)

def partitionElements(Comm comm, int nLayersOfOverlap, cmeshToolsTest.CMesh cmesh, cmeshToolsTest.CMesh subdomain_cmesh):
    cmesh.meshlink.mesh.subdomainp = &subdomain_cmesh.meshlink.mesh
    c_partitionElements(comm.ob_mpi,
                        cmesh.meshlink.mesh,
                        nLayersOfOverlap)
    print("infwt nNodes_global={0:d}".format(cmesh.meshlink.mesh.nNodes_global))
    return (
        np.asarray(<int[:(comm.size+1)]> cmesh.meshlink.mesh.elementOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.meshlink.mesh.subdomainp.nElements_global]> cmesh.meshlink.mesh.elementNumbering_subdomain2global),
        np.asarray(<int[:(comm.size+1)]> cmesh.meshlink.mesh.nodeOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.meshlink.mesh.subdomainp.nNodes_global]> cmesh.meshlink.mesh.nodeNumbering_subdomain2global),
        np.asarray(<int[:(comm.size+1)]> cmesh.meshlink.mesh.elementBoundaryOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.meshlink.mesh.subdomainp.nElementBoundaries_global]> cmesh.meshlink.mesh.elementBoundaryNumbering_subdomain2global),
        np.asarray(<int[:(comm.size+1)]> cmesh.meshlink.mesh.edgeOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.meshlink.mesh.subdomainp.nEdges_global]> cmesh.meshlink.mesh.edgeNumbering_subdomain2global)
    )

def partitionNodes(Comm comm, int nLayersOfOverlap, cmeshToolsTest.CMesh cmesh, cmeshToolsTest.CMesh subdomain_cmesh):
    cmesh.meshlink.mesh.subdomainp = &subdomain_cmesh.meshlink.mesh
    c_partitionNodes(comm.ob_mpi,
                     cmesh.meshlink.mesh,
                     nLayersOfOverlap)
    print("infwt nNodes_global={0:d}".format(cmesh.meshlink.mesh.nNodes_global))
    return (
        np.asarray(<int[:comm.size+1]> cmesh.meshlink.mesh.elementOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.meshlink.mesh.subdomainp.nElements_global]> cmesh.meshlink.mesh.elementNumbering_subdomain2global),
        np.asarray(<int[:comm.size+1]> cmesh.meshlink.mesh.nodeOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.meshlink.mesh.subdomainp.nNodes_global]> cmesh.meshlink.mesh.nodeNumbering_subdomain2global),
        np.asarray(<int[:comm.size+1]> cmesh.meshlink.mesh.elementBoundaryOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.meshlink.mesh.subdomainp.nElementBoundaries_global]> cmesh.meshlink.mesh.elementBoundaryNumbering_subdomain2global),
        np.asarray(<int[:comm.size+1]> cmesh.meshlink.mesh.edgeOffsets_subdomain_owned),
        np.asarray(<int[:cmesh.meshlink.mesh.subdomainp.nEdges_global]> cmesh.meshlink.mesh.edgeNumbering_subdomain2global)
    )

