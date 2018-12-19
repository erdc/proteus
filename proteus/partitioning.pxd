# A type of -*- python -*- file
from mpi4py.libmpi cimport MPI_Comm
cimport mesh

cdef extern from "partitioning.h" namespace "proteus":
    extern int c_partitionNodes "proteus::partitionNodes" (const MPI_Comm& PROTEUS_COMM_WORLD,
                                                           mesh.Mesh& mesh,
                                                           int nNodes_overlap);
    extern int c_partitionNodesFromTetgenFiles "proteus::partitionNodesFromTetgenFiles" (const MPI_Comm& PROTEUS_COMM_WORLD,
                                                                                         const char* filebase,
                                                                                         int indexBase, mesh.Mesh& newMesh,
                                                                                         int nNodes_overlap);
    extern int c_partitionElements "proteus::partitionElements" (const MPI_Comm& PROTEUS_COMM_WORLD,
                                                                 mesh.Mesh& mesh,
                                                                 int nElements_overlap);
    extern int buildQuadraticSubdomain2GlobalMappings_1d(const MPI_Comm& PROTEUS_COMM_WORLD, mesh.Mesh& mesh,
                                                         const int *elementOffsets_subdomain_owned,
                                                         const int *nodeOffsets_subdomain_owned,
                                                         const int *elementNumbering_subdomain2global,
                                                         const int *nodeNumbering_subdomain2global,
                                                         int& nDOF_all_processes,#total number of dofs in whole domain
                                                         int& nDOF_subdomain,#total number of dofs in sub-domain
                                                         int& max_dof_neighbors,#maximum number of neighbors for connectivity of dofs
                                                         int *offsets_subdomain_owned, #starting point of local dofs on each processor (nProcs+1)
                                                         int * subdomain_l2g, #local to global dof mapping on subdomain
                                                         int* subdomain2global,#subdomain dof to global (parallel) numbering
                                                         double * lagrangeNodesArray);#location of nodes corresponding to dofs
    extern int buildQuadraticSubdomain2GlobalMappings_2d(const MPI_Comm& PROTEUS_COMM_WORLD, mesh.Mesh& mesh,
                                                         const int *elementBoundaryOffsets_subdomain_owned,
                                                         const int *nodeOffsets_subdomain_owned,
                                                         const int *elementBoundaryNumbering_subdomain2global,
                                                         const int *nodeNumbering_subdomain2global,
                                                         int& nDOF_all_processes,#total number of dofs in whole domain
                                                         int& nDOF_subdomain,#total number of dofs in sub-domain
                                                         int& max_dof_neighbors,#maximum number of neighbors for connectivity of dofs
                                                         int *offsets_subdomain_owned, #starting point of local dofs on each processor (nProcs+1)
                                                         int *subdomain_l2g, #local to global dof mapping on subdomain
                                                         int *subdomain2global,#subdomain dof to global (parallel) numbering
                                                         double * lagrangeNodesArray);#location of nodes corresponding to dofs
    extern int buildQuadraticSubdomain2GlobalMappings_3d(const MPI_Comm& PROTEUS_COMM_WORLD, mesh.Mesh& mesh,
                                                         const int *edgeOffsets_subdomain_owned,
                                                         const int *nodeOffsets_subdomain_owned,
                                                         const int *edgeNumbering_subdomain2global,
                                                         const int *nodeNumbering_subdomain2global,
                                                         int& nDOF_all_processes,#total number of dofs in whole domain
                                                         int& nDOF_subdomain,#total number of dofs in sub-domain
                                                         int& max_dof_neighbors,#maximum number of neighbors for connectivity of dofs
                                                         int *offsets_subdomain_owned, #starting point of local dofs on each processor (nProcs+1)
                                                         int *subdomain_l2g, #local to global dof mapping on subdomain
                                                         int *subdomain2global,#subdomain dof to global (parallel) numbering
                                                         double * lagrangeNodesArray);#location of nodes corresponding to dofs
    extern int buildQuadraticCubeSubdomain2GlobalMappings_3d(const MPI_Comm& PROTEUS_COMM_WORLD, mesh.Mesh& mesh,
                                                             const int *edgeOffsets_subdomain_owned,
                                                             const int *nodeOffsets_subdomain_owned,
                                                             const int *edgeNumbering_subdomain2global,
                                                             const int *nodeNumbering_subdomain2global,
                                                             int& nDOF_all_processes,#total number of dofs in whole domain
                                                             int& nDOF_subdomain,#total number of dofs in sub-domain
                                                             int& max_dof_neighbors,#maximum number of neighbors for connectivity of dofs
                                                             int *offsets_subdomain_owned, #starting point of local dofs on each processor (nProcs+1)
                                                             int *subdomain_l2g, #local to global dof mapping on subdomain
                                                             int *subdomain2global,#subdomain dof to global (parallel) numbering
                                                             double * lagrangeNodesArray);#location of nodes corresponding to dofs
    extern int buildDiscontinuousGalerkinSubdomain2GlobalMappings(const MPI_Comm& PROTEUS_COMM_WORLD, mesh.Mesh& mesh,
                                                                  const int *elementOffsets_subdomain_owned,
                                                                  const int *elementNumbering_subdomain2global,
                                                                  int nDOF_element,
                                                                  int& nDOF_all_processes,#total number of dofs in whole domain
                                                                  int& nDOF_subdomain,#total number of dofs in sub-domain
                                                                  int& max_dof_neighbors,#maximum number of neighbors for connectivity of dofs
                                                                  int *offsets_subdomain_owned, #starting point of local dofs on each processor (nProcs+1)
                                                                  int * subdomain_l2g, #local to global dof mapping on subdomain
                                                                  int* subdomain2global);#subdomain dof to global (parallel) numbering
