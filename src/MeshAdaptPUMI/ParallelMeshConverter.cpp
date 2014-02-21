#include <algorithm>

#include "MeshAdaptPUMI.h"
#include "mesh.h"

const int DEFAULT_ELEMENT_MATERIAL=0;
const int DEFAULT_NODE_MATERIAL=-1;
const int INTERIOR_NODE_MATERIAL=0;
const int EXTERIOR_NODE_MATERIAL=1;
const int INTERIOR_ELEMENT_BOUNDARY_MATERIAL=0;
const int EXTERIOR_ELEMENT_BOUNDARY_MATERIAL=1;

int MeshAdaptPUMIDrvr::ConstructFromParallelPUMIMesh(Mesh& mesh, Mesh& subdomain_mesh) {

//Current impementation is for serial only. Need to resolve things with Proteus developers for parallel implementation  
  
  mesh.subdomainp = new Mesh();

  mesh.subdomainp = &subdomain_mesh;
  std::cout << "Constructing global data structures\n"; 
  int numGlobElem;
  PUMI_Mesh_GetNumEnt(PUMI_MeshInstance, PUMI_REGION, PUMI_ALLTOPO, &numGlobElem);
  mesh.nElements_global = numGlobElem;

  int numLocElem;
  PUMI_Part_GetNumEnt(PUMI_Part, PUMI_REGION, PUMI_ALLTOPO, &numLocElem);
  mesh.subdomainp->nElements_global = numLocElem;

  int numGlobNodes;
  PUMI_Mesh_GetNumEnt(PUMI_MeshInstance, PUMI_VERTEX, PUMI_ALLTOPO, &numGlobNodes);
  mesh.nNodes_global = numGlobNodes;

  int numLocNodes;
  PUMI_Part_GetNumEnt(PUMI_Part, PUMI_VERTEX, PUMI_ALLTOPO, &numLocNodes);
  mesh.subdomainp->nNodes_global = numLocNodes;

  int numGlobFaces;
  PUMI_Mesh_GetNumEnt(PUMI_MeshInstance, PUMI_FACE, PUMI_ALLTOPO, &numGlobFaces);
  mesh.nElementBoundaries_global = numGlobFaces;

  int numLocFaces;
  PUMI_Part_GetNumEnt(PUMI_Part, PUMI_FACE, PUMI_ALLTOPO, &numLocFaces);
  mesh.subdomainp->nElementBoundaries_global = numLocFaces;

  int numGlobEdges;
  PUMI_Mesh_GetNumEnt(PUMI_MeshInstance, PUMI_EDGE, PUMI_ALLTOPO, &numGlobEdges);
  mesh.nEdges_global = numGlobEdges;

  int numLocEdges;
  PUMI_Part_GetNumEnt(PUMI_Part, PUMI_EDGE, PUMI_ALLTOPO, &numLocEdges);
  mesh.subdomainp->nEdges_global = numLocEdges;
//nNodes_element for now is constant for the entire mesh, Ask proteus about using mixed meshes
//therefore currently this code only supports tet meshes
  mesh.subdomainp->nNodes_element = 4; //hardcore: for tets, number of nodes per element
  mesh.subdomainp->nNodes_elementBoundary = 3; //hardcode: for tets, looks like number of nodes of a face
  mesh.subdomainp->nElementBoundaries_element = 4; //hardcode: for tets, looks like number of faces
  
  std::cerr << "*******Local Proteus Mesh Stats*********\n";
  std::cerr << "Rank: " << comm_rank << ": Number of elements " << mesh.subdomainp->nElements_global << "\n";
  std::cerr << "Rank: " << comm_rank << ": Number of nodes " << mesh.subdomainp->nNodes_global << "\n";
  std::cerr << "Rank: " << comm_rank << ": Number of boundaries " << mesh.subdomainp->nElementBoundaries_global << "\n";
  std::cerr << "Rank: " << comm_rank << ": Number of edges " << mesh.subdomainp->nEdges_global << "\n";
  std::cerr << "*****************************************\n";
  
  if(comm_rank==0) {
    std::cerr << "*******Global Proteus Mesh Stats*********\n";
    std::cerr << ": Number of elements " << mesh.nElements_global << "\n";
    std::cerr << ": Number of nodes " << mesh.nNodes_global << "\n";
    std::cerr << ": Number of boundaries " << mesh.nElementBoundaries_global << "\n";
    std::cerr << ": Number of edges " << mesh.nEdges_global << "\n";
    std::cerr << "*****************************************\n";
  }

  ConstructGlobalNumbering(mesh);
  ConstructNodes(*mesh.subdomainp);
  ConstructElements(*mesh.subdomainp);
  ConstructBoundaries(*mesh.subdomainp);
  ConstructEdges(*mesh.subdomainp);
  ConstructGlobalStructures(mesh); 

  ConstructMaterialArrays(*mesh.subdomainp);

   //chitak: debug
/*  
//debugging and checking
   for(int i=0; i<mesh.subdomainp->nElements_global; i++) {
     std::cout << "element, local: " << i << " global: " << mesh.elementNumbering_subdomain2global[i] << "\n";
   }
   MPI_Barrier(MPI_COMM_WORLD);
   for(int i=0; i<mesh.subdomainp->nNodes_global; i++) {
     std::cout << "nodes, local: " << i << " global: " << mesh.nodeNumbering_subdomain2global[i] << "\n";
   }
   for(int i=0; i<comm_size+1; i++) {
      std::cout << "my rank: " << comm_rank << " element_subdomain_owned: " << mesh.elementOffsets_subdomain_owned[i] << "\n";
   }
   for(int i=0; i<comm_size+1; i++) {
      std::cout << "my rank: " << comm_rank << " node_subdomain_owned: " << mesh.nodeOffsets_subdomain_owned[i] << "\n";
   }
*/
  return 0;
} 

int MeshAdaptPUMIDrvr::ConstructGlobalStructures(Mesh &mesh) {

  mesh.elementNumbering_subdomain2global = new int[mesh.subdomainp->nElements_global];
  mesh.elementBoundaryNumbering_subdomain2global = new int[mesh.subdomainp->nElementBoundaries_global];
  mesh.nodeNumbering_subdomain2global = new int[mesh.subdomainp->nNodes_global];
  mesh.edgeNumbering_subdomain2global = new int[mesh.subdomainp->nEdges_global];

  int partDim;
  PUMI_Part_GetDim(PUMI_Part, &partDim);
  for (int type=partDim; type>-1; type--) {
    
    int* temp_subdomain2global;
    if(type==3) temp_subdomain2global = mesh.elementNumbering_subdomain2global;
    if(type==2) temp_subdomain2global = mesh.elementBoundaryNumbering_subdomain2global;
    if(type==1) temp_subdomain2global = mesh.edgeNumbering_subdomain2global;
    if(type==1) temp_subdomain2global = mesh.nodeNumbering_subdomain2global;
    
    int isEnd = 0;
    int eN = 0;
    pPartEntIter EntIt;
    pMeshEnt meshEnt; 
    PUMI_PartEntIter_Init (PUMI_Part, type, PUMI_ALLTOPO, EntIt);
    while (!isEnd) {
      PUMI_PartEntIter_GetNext(EntIt, meshEnt);

      int entID = PUMI_MeshEnt_ID(meshEnt);
      int globID;
      PUMI_MeshEnt_GetIntTag(PUMI_MeshInstance, meshEnt, GlobNumberTag, &globID);
      
      temp_subdomain2global[entID] = globID;

      PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
    } //region loop
    PUMI_PartEntIter_Del(EntIt);
  } //type


}

int MeshAdaptPUMIDrvr::ConstructGlobalNumbering(Mesh &mesh) {
 
 PUMI_Mesh_CreateTag(PUMI_MeshInstance, "GlobNumbering", PUMI_INT, 1, GlobNumberTag);
 
 mesh.elementOffsets_subdomain_owned = new int[comm_size+1];
 mesh.elementBoundaryOffsets_subdomain_owned = new int[comm_size+1];
 mesh.edgeOffsets_subdomain_owned = new int[comm_size+1];
 mesh.nodeOffsets_subdomain_owned = new int[comm_size+1];

 int partDim;
 PUMI_Part_GetDim(PUMI_Part, &partDim);
 for (int type=partDim; type>-1; type--) {
  
   int nLocalOwned;
   if(type==3){ //elements, we own all elements on a part
     PUMI_Part_GetNumEnt(PUMI_Part, type, PUMI_ALLTOPO, &nLocalOwned); //we don't have overlapping elements, all elements are owned on a part
   } else { //for vertices, edges, faces
     CalculateOwnedEnts(static_cast<PUMI_EntType>(type), nLocalOwned);
   }

   int nOwned[comm_size];
   CommunicateOwnedNumbers(nLocalOwned, nOwned);

   int* numOffset;
   if(type==3){
     numOffset = mesh.elementOffsets_subdomain_owned;
     elms_owned = nLocalOwned;
   }
   if(type==2){
     numOffset = mesh.elementBoundaryOffsets_subdomain_owned;
     faces_owned = nLocalOwned;
   }
   if(type==1){
     numOffset = mesh.edgeOffsets_subdomain_owned;
     edges_owned = nLocalOwned;
   }
   if(type==0){
     numOffset = mesh.nodeOffsets_subdomain_owned;
     vtx_owned = nLocalOwned;
   }

   numOffset[0]=0;
   for(int i=1; i<comm_size+1; i++) {
     numOffset[i] = numOffset[i-1]+nOwned[i-1];
   }
   for(int i=0; i<comm_size; i++) {
      std::cout << "my rank: " << comm_rank <<" rank: " << i << " : Entity " << type << ": " << nOwned[i] << "\n";
   }
 
   SetOwnerGlobNumbering(GlobNumberTag, static_cast<PUMI_EntType>(type), numOffset[comm_rank]);
   if(type!=3) 
     SetCopyGlobNumbering(GlobNumberTag, type);

   //communicate numOffset to populate mesh.offset_owned etc

 } //loop on entity types

/* 
//debugging and checking :
  pPartEntIter EntIt;
  pMeshEnt meshEnt;
  PUMI_PartEntIter_Init (PUMI_Part, PUMI_VERTEX, PUMI_ALLTOPO, EntIt);

  int isEnd = 0;
  int eN = 0;
  while (!isEnd) {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);

    int globNum;
    PUMI_MeshEnt_GetIntTag(PUMI_MeshInstance, meshEnt, GlobNumberTag, &globNum);
    std::cout << "rank: " << comm_rank << ": local num: " << eN << " Glob num: " << globNum << "\n";
    PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
    eN++;
  }

  PUMI_PartEntIter_Del(EntIt);
*/
 return 0; 
}

int MeshAdaptPUMIDrvr::CalculateOwnedEnts(PUMI_EntType EntType, int &nOwnedEnt ) {  

  pPartEntIter EntIt;
  pMeshEnt meshEnt;
  PUMI_PartEntIter_Init (PUMI_Part, EntType, PUMI_ALLTOPO, EntIt);

  int isEnd = 0;
  int eN = 0;
  
  int nLocalOwnedNodes;
  while (!isEnd) {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);

    int owned;
    PUMI_MeshEnt_IsOwned(meshEnt, PUMI_Part, &owned);
    if(owned) {
       eN++;
    }

    PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
  } 
   
  PUMI_PartEntIter_Del(EntIt);
  nOwnedEnt = eN;

  return 0;
}

int MeshAdaptPUMIDrvr::CommunicateOwnedNumbers(int toSend, int *toReceive) {
  
  MPI_Status status;

  MPI_Send(&toSend, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); //send local number of owned to process 0
  
  if(comm_rank == 0) {
    toReceive[0] = toSend;
    for(int i=1; i<comm_size; i++) {
      MPI_Recv(&toReceive[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
    }
  }   

  MPI_Bcast(&toReceive[0], comm_size, MPI_INT, 0, MPI_COMM_WORLD);
  
  return 0;
}

int MeshAdaptPUMIDrvr::SetOwnerGlobNumbering(pTag globNumTag, PUMI_EntType EntType, int numOffset) {
  
  pPartEntIter EntIt;
  pMeshEnt meshEnt;
  PUMI_PartEntIter_Init(PUMI_Part, EntType, PUMI_ALLTOPO, EntIt);

  int isEnd = 0;
  int eN = 0;
  while (!isEnd) {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);

    int owned;
    PUMI_MeshEnt_IsOwned(meshEnt, PUMI_Part, &owned);
    if(owned) {
//       locGlobNumbering[numOffset+eN] = eN; //global element numbering 
       PUMI_MeshEnt_SetIntTag(PUMI_MeshInstance, meshEnt, globNumTag, numOffset+eN);
//       PUMI_MeshEnt_SetID(meshEnt, numOffset+eN);
       eN++;
    }

    PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
  } //region loop

  PUMI_PartEntIter_Del(EntIt);
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}

int MeshAdaptPUMIDrvr::SetCopyGlobNumbering(pTag globNumTag, int EntType) {

 if (SCUTIL_CommSize() > 1) { //need to do this only for parallel

  pPartEntIter EntIt;
  pMeshEnt meshEnt;
  PUMI_PartEntIter_InitPartBdry(PUMI_Part, -1, EntType, PUMI_ALLTOPO, EntIt);
  int isEnd = 0;
  int eN = 0;

  PCU_Comm_Start (PCU_GLOBAL_METHOD);
  int num_sent = 0, num_recvd = 0;
//  std::pair<pMeshEnt, int>* msg_send = (std::pair<pMeshEnt, int>)malloc(sizeof(std::pair<pMeshEnt, int>));
  std::pair<pMeshEnt, int>* msg_send = new std::pair<pMeshEnt, int>;
  size_t msg_size = sizeof(std::pair<pMeshEnt,int>);

  while (!isEnd) {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);

    int owned;
    PUMI_MeshEnt_IsOwned(meshEnt, PUMI_Part, &owned);
    if(owned) {
      int ownGlobNum;
      PUMI_MeshEnt_GetIntTag(PUMI_MeshInstance, meshEnt, globNumTag, &ownGlobNum);
      msg_send->second =  ownGlobNum;

      std::vector<std::pair<int, pMeshEnt> > VecRemCpy;
      PUMI_MeshEnt_GetAllRmt(meshEnt, VecRemCpy);
      for (int iRC = 0; iRC < VecRemCpy.size(); ++iRC)
      {
        msg_send->first = VecRemCpy[iRC].second;
        PCU_Comm_Write (VecRemCpy[iRC].first, (void*)msg_send, msg_size);
      }
    }
    PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
    eN++;
  }
  PUMI_PartEntIter_Del(EntIt);
  MPI_Barrier(MPI_COMM_WORLD);

  delete msg_send;
  PCU_Comm_Send ();

  size_t recv_size;
  int pid_from;
  void* msg_recv;
  while (PCU_Comm_Read (&pid_from,&msg_recv,&recv_size)) {
     std::pair<pMeshEnt, int>& msg_pair = *(static_cast<std::pair<pMeshEnt, int> *>(msg_recv));
     pMeshEnt copyEnt = msg_pair.first;
     PUMI_MeshEnt_SetIntTag(PUMI_MeshInstance, copyEnt, globNumTag, msg_pair.second);
  }

 } //if parallel
  return 0;

}
