#include <algorithm>
#include <mpi.h>
#include <PCU.h>
#include <sstream>

#include "MeshAdaptPUMI.h"
#include "mesh.h"

#include "DumpMesh.h"

const int DEFAULT_ELEMENT_MATERIAL=0;
const int DEFAULT_NODE_MATERIAL=-1;
const int INTERIOR_NODE_MATERIAL=0;
const int EXTERIOR_NODE_MATERIAL=1;
const int INTERIOR_ELEMENT_BOUNDARY_MATERIAL=0;
const int EXTERIOR_ELEMENT_BOUNDARY_MATERIAL=1;

static int countTotal(apf::Mesh* m, int dim)
{
  int total = apf::countOwned(m, dim);
  PCU_Add_Ints(&total, 1);
  return total;
}

int MeshAdaptPUMIDrvr::constructFromParallelPUMIMesh(Mesh& mesh, Mesh& subdomain_mesh)
{
  mesh.subdomainp = new Mesh();
  mesh.subdomainp = &subdomain_mesh;
  initializeMesh(subdomain_mesh);
  if (!PCU_Comm_Self())
    std::cerr << "Constructing parallel proteus mesh\n"; 
  
  int numGlobElem = countTotal(m, 3);
  mesh.nElements_global = numGlobElem;

  int numLocElem = m->count(3);
  mesh.subdomainp->nElements_global = numLocElem;

  int numGlobNodes = countTotal(m, 0);
  mesh.nNodes_global = numGlobNodes;

  int numLocNodes = m->count(0);
  mesh.subdomainp->nNodes_global = numLocNodes;

  int numGlobFaces = countTotal(m, 2);
  mesh.nElementBoundaries_global = numGlobFaces;

  int numLocFaces = m->count(2);
  mesh.subdomainp->nElementBoundaries_global = numLocFaces;

  int numGlobEdges = countTotal(m, 1);
  mesh.nEdges_global = numGlobEdges;

  int numLocEdges = m->count(1);
  mesh.subdomainp->nEdges_global = numLocEdges;
//nNodes_element for now is constant for the entire mesh, Ask proteus about using mixed meshes
//therefore currently this code only supports tet meshes
  mesh.subdomainp->nNodes_element = 4; //hardcode: for tets, number of nodes per element
  mesh.subdomainp->nNodes_elementBoundary = 3; //hardcode: for tets, looks like number of nodes of a face
  mesh.subdomainp->nElementBoundaries_element = 4; //hardcode: for tets, looks like number of faces/element
  
#ifdef MESH_INFO
  for (int i = 0; i < PCU_Comm_Peers(); ++i) {
    if (i == PCU_Comm_Self()) {
      std::cerr << "*******Local Proteus Mesh Stats*********\n";
      std::cerr << "Rank: " << comm_rank << ": Number of elements " << mesh.subdomainp->nElements_global << "\n";
      std::cerr << "Rank: " << comm_rank << ": Number of nodes " << mesh.subdomainp->nNodes_global << "\n";
      std::cerr << "Rank: " << comm_rank << ": Number of boundaries " << mesh.subdomainp->nElementBoundaries_global << "\n";
      std::cerr << "Rank: " << comm_rank << ": Number of edges " << mesh.subdomainp->nEdges_global << "\n";
      std::cerr << "*****************************************\n";
    }
    PCU_Barrier();
  }
  
  if(comm_rank==0) {
    std::cerr << "*******Global Proteus Mesh Stats*********\n";
    std::cerr << ": Number of elements " << mesh.nElements_global << "\n";
    std::cerr << ": Number of nodes " << mesh.nNodes_global << "\n";
    std::cerr << ": Number of boundaries " << mesh.nElementBoundaries_global << "\n";
    std::cerr << ": Number of edges " << mesh.nEdges_global << "\n";
    std::cerr << "*****************************************\n";
  }
#endif

  constructGlobalNumbering(mesh);
  numberLocally();
  constructNodes(*mesh.subdomainp);
  constructElements(*mesh.subdomainp);
  constructBoundaries(*mesh.subdomainp);
  constructEdges(*mesh.subdomainp);
  constructMaterialArrays(*mesh.subdomainp);
  constructGlobalStructures(mesh); 

  std::stringstream ss;
  ss << "ToProteus_t" << nAdapt<<".smb";
  std::string s = ss.str();
  m->writeNative(s.c_str());

  return 0;
} 

int MeshAdaptPUMIDrvr::constructGlobalNumbering(Mesh &mesh)
{
  /* N^2 data structures and algorithms are terrible for scalability.
     Going along because Proteus is structured this way. */
  mesh.elementOffsets_subdomain_owned = new int[comm_size+1];
  mesh.elementBoundaryOffsets_subdomain_owned = new int[comm_size+1];
  mesh.edgeOffsets_subdomain_owned = new int[comm_size+1];
  mesh.nodeOffsets_subdomain_owned = new int[comm_size+1];

  for (int dim = 0; dim <= m->getDimension(); ++dim) {

    int nLocalOwned = apf::countOwned(m, dim);
    int localOffset = nLocalOwned;
    PCU_Exscan_Ints(&localOffset, 1);

    int* allOffsets;
    if(dim==3){
      allOffsets = mesh.elementOffsets_subdomain_owned;
    }
    if(dim==2){
      allOffsets = mesh.elementBoundaryOffsets_subdomain_owned;
    }
    if(dim==1){
      allOffsets = mesh.edgeOffsets_subdomain_owned;
    }
    if(dim==0){
      allOffsets = mesh.nodeOffsets_subdomain_owned;
    }

    /* one of the many reasons N^2 algorithms are bad */
    MPI_Allgather(&localOffset, 1, MPI_INT,
                  allOffsets, 1, MPI_INT, MPI_COMM_WORLD);
    allOffsets[PCU_Comm_Peers()] = countTotal(m, dim);

    std::stringstream ss;
    ss << "proteus_global_";
    ss << dim;
    std::string name = ss.str();
    /* this algorithm does global numbering properly,
       without O(#procs) runtime */
    global[dim] = apf::makeGlobal(apf::numberOwnedDimension(m, name.c_str(), dim));
    apf::synchronize(global[dim]);
  } //loop on entity dimensions

  return 0; 
}

int MeshAdaptPUMIDrvr::constructGlobalStructures(Mesh &mesh)
{
  mesh.elementNumbering_subdomain2global = new int[mesh.subdomainp->nElements_global];
  mesh.elementBoundaryNumbering_subdomain2global = new int[mesh.subdomainp->nElementBoundaries_global];
  mesh.nodeNumbering_subdomain2global = new int[mesh.subdomainp->nNodes_global];
  mesh.edgeNumbering_subdomain2global = new int[mesh.subdomainp->nEdges_global];
  
  for (int d = 0; d <= m->getDimension(); ++d) {
    int* temp_subdomain2global;
    if(d==3) temp_subdomain2global = mesh.elementNumbering_subdomain2global;
    if(d==2) temp_subdomain2global = mesh.elementBoundaryNumbering_subdomain2global;
    if(d==1) temp_subdomain2global = mesh.edgeNumbering_subdomain2global;
    if(d==0) temp_subdomain2global = mesh.nodeNumbering_subdomain2global;

    apf::MeshIterator* it = m->begin(d);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
      int i = localNumber(e);
      temp_subdomain2global[i] = apf::getNumber(global[d], apf::Node(e, 0));
    }
    m->end(it);
    apf::destroyGlobalNumbering(global[d]);
  }

  return 0;
}

int MeshAdaptPUMIDrvr::dumpMesh(Mesh& mesh)
{
  std::stringstream ss;
  ss << "dan_debug_" << PCU_Comm_Self() << ".txt";
  std::string s = ss.str();
  FILE* f = fopen(s.c_str(), "w");
  dump_proteus_mesh(&mesh, f);
  fclose(f);
  return 0;
}

