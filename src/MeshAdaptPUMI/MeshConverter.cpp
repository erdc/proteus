#include <algorithm>

#include "MeshAdaptPUMI.h"
#include "mesh.h"

const int DEFAULT_ELEMENT_MATERIAL=0;
const int DEFAULT_NODE_MATERIAL=-1;
const int INTERIOR_NODE_MATERIAL=0;
const int EXTERIOR_NODE_MATERIAL=1;
const int INTERIOR_ELEMENT_BOUNDARY_MATERIAL=0;
const int EXTERIOR_ELEMENT_BOUNDARY_MATERIAL=1;

int MeshAdaptPUMIDrvr::ConstructFromSerialPUMIMesh(Mesh& mesh) {

//  Mesh mesh = MESH(cmesh);

//Current impementation is for serial only. Need to resolve things with Proteus developers for parallel implementation  
  std::cout << "Constructing global data structures\n"; 
  int numGlobElem;
  PUMI_Mesh_GetNumEnt(PUMI_MeshInstance, PUMI_REGION, PUMI_ALLTOPO, &numGlobElem);
  mesh.nElements_global = numGlobElem;
  elms_owned = numGlobElem;

  int numGlobNodes;
  PUMI_Mesh_GetNumEnt(PUMI_MeshInstance, PUMI_VERTEX, PUMI_ALLTOPO, &numGlobNodes);
  mesh.nNodes_global = numGlobNodes;
  vtx_owned = numGlobNodes;

  int numGlobFaces;
  PUMI_Mesh_GetNumEnt(PUMI_MeshInstance, PUMI_FACE, PUMI_ALLTOPO, &numGlobFaces);
  mesh.nElementBoundaries_global = numGlobFaces;
  faces_owned = numGlobFaces;

  int numGlobEdges;
  PUMI_Mesh_GetNumEnt(PUMI_MeshInstance, PUMI_EDGE, PUMI_ALLTOPO, &numGlobEdges);
  mesh.nEdges_global = numGlobEdges;
  edges_owned = numGlobEdges;

//nNodes_element for now is constant for the entire mesh, Ask proteus about using mixed meshes
//therefore currently this code only supports tet meshes
  mesh.nNodes_element = 4; //hardcore: for tets, number of nodes per element
  mesh.nNodes_elementBoundary = 3; //hardcode: for tets, looks like number of nodes of a face
  mesh.nElementBoundaries_element = 4; //hardcode: for tets, looks like number of faces

  std::cerr << "*******Proteus Mesh Stats*********\n";
  std::cerr << "Number of elements " << mesh.nElements_global << "\n";
  std::cerr << "Number of nodes " << mesh.nNodes_global << "\n";
  std::cerr << "Number of boundaries " << mesh.nElementBoundaries_global << "\n";
  std::cerr << "Number of edges " << mesh.nEdges_global << "\n";
  
  ConstructNodes(mesh);
  ConstructElements(mesh);
  ConstructBoundaries(mesh);
  ConstructEdges(mesh);
  ConstructMaterialArrays(mesh);

  return 0;
}  

int MeshAdaptPUMIDrvr::ConstructNodes(Mesh& mesh) {
//////build node array (coordinates?)

  std::cout << "Constructing nodal data structures\n" ; 
//set interior and exterior node material flags after get boundary info in
//constructElementBoundaryElementsArray_edge 
//if nodeMaterialTypes left as DEFAULT_NODE_MATERIALi
  mesh.nodeArray = new double[mesh.nNodes_global*3];
  memset(mesh.nodeArray,0,mesh.nNodes_global*3*sizeof(double));
  mesh.nodeMaterialTypes = new int[mesh.nNodes_global];
  memset(mesh.nodeMaterialTypes,DEFAULT_NODE_MATERIAL,mesh.nNodes_global*sizeof(int));

  pPartEntIter EntIt;
  pMeshEnt meshEnt;
  PUMI_PartEntIter_Init (PUMI_Part, PUMI_VERTEX, PUMI_ALLTOPO, EntIt);
  int isEnd = 0;
  int nN = 0; int nN_owned = 0; 
  int nN_copy = vtx_owned;
  while (!isEnd) {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);
    
    int owned;
    PUMI_MeshEnt_IsOwned(meshEnt, PUMI_Part, &owned);
    if(owned) {
       PUMI_MeshEnt_SetID(meshEnt, nN_owned);
       nN_owned++;
    } else {
       PUMI_MeshEnt_SetID(meshEnt, nN_copy);
       nN_copy++;
    }
    
    int VtxID = PUMI_MeshEnt_ID(meshEnt);
    double coords[3]={0.,0.,0.};
    double *xyz = coords;
    PUMI_MeshVtx_GetCoord(meshEnt, xyz);
    
    for(int i=0; i<3; i++) {
      mesh.nodeArray[VtxID*3+i]= coords[i];
    }
    nN++;
    PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
  }
  PUMI_PartEntIter_Del(EntIt);
  assert(nN_owned==vtx_owned-1);
  assert(nN_copy==mesh.nNodes_global-1);
  
  return 0;
//
//end: build node array
//
} 

int MeshAdaptPUMIDrvr::ConstructElements(Mesh& mesh) {

  std::cout << "Constructing element data structures\n"; 
//////build node list for the elements
  mesh.elementNodesArray = new int[mesh.nElements_global*mesh.nNodes_element];
  mesh.elementMaterialTypes = new int[mesh.nElements_global];
  memset(mesh.elementMaterialTypes,DEFAULT_ELEMENT_MATERIAL,mesh.nElements_global*sizeof(int));
//
  pPartEntIter EntIt;
  pMeshEnt meshEnt;
  PUMI_PartEntIter_Init (PUMI_Part, PUMI_REGION, PUMI_ALLTOPO, EntIt);

  int isEnd = 0;
  int eN = 0; int eN_owned = 0;
  int eN_copy = elms_owned;
  while (!isEnd) {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);

    int owned;
    PUMI_MeshEnt_IsOwned(meshEnt, PUMI_Part, &owned);
    if(owned) {
       PUMI_MeshEnt_SetID(meshEnt, eN_owned);
       eN_owned++;
    } else {
       PUMI_MeshEnt_SetID(meshEnt, eN_copy);
       eN_copy++;
    }
    
    int elemID = PUMI_MeshEnt_ID(meshEnt);
    std::vector<pMeshEnt> vecVtx;
    PUMI_MeshEnt_GetAdj(meshEnt, PUMI_VERTEX, 0, vecVtx);
    int iNumVtx = vecVtx.size();
    for (int iVtx = 0; iVtx < iNumVtx; ++iVtx) {
      int VtxID = PUMI_MeshEnt_ID(vecVtx[iVtx]);
      mesh.elementNodesArray[elemID*mesh.nNodes_element+iVtx] = VtxID;  //populating the element nodes array 
//      std::cerr << "Element nodes array: " << mesh.elementNodesArray[eN*mesh.nNodes_element+iVtx]<< "\n";
    }
    PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
    eN++; //increment element number of the array
  } //region loop

  PUMI_PartEntIter_Del(EntIt);
  assert(eN_owned==elms_owned-1);
  assert(eN_copy==mesh.nElements_global-1);

  return 0;
//end: build node list for the elements
}

int MeshAdaptPUMIDrvr::ConstructBoundaries(Mesh& mesh) {

  std::cout << "Constructing boundary data structures\n" ; 
////////build face list (elementBoundary and nodeBoundary arrays) : this is a bit complex and should be verified
//Enter at your own peril for those who stray will be lost
//
  std::set<int> interiorElementBoundaries,exteriorElementBoundaries; //do we need this? we can do this by looping over model faces, right?

  mesh.elementBoundaryNodesArray =  new int[mesh.nElementBoundaries_global*mesh.nNodes_elementBoundary];
  mesh.elementBoundaryElementsArray = new int[mesh.nElementBoundaries_global*2];
  mesh.elementBoundaryLocalElementBoundariesArray = new int[mesh.nElementBoundaries_global*2];
  mesh.elementNeighborsArray = new int[mesh.nElements_global*mesh.nElementBoundaries_element];
  //mwf added
  mesh.elementBoundariesArray= new int[mesh.nElements_global*mesh.nElementBoundaries_element];

  pPartEntIter EntIt;
  pMeshEnt meshEnt;
  PUMI_PartEntIter_Init (PUMI_Part, PUMI_FACE, PUMI_TRI, EntIt); //currently only implemented for tetrahedra
  int isEnd = 0;
  int nF = 0; int nF_owned = 0;
  int nF_copy = faces_owned;
  while (!isEnd) {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);

    int owned;
    PUMI_MeshEnt_IsOwned(meshEnt, PUMI_Part, &owned);
    if(owned) {
       PUMI_MeshEnt_SetID(meshEnt, nF_owned);
       nF_owned++;
    } else {
       PUMI_MeshEnt_SetID(meshEnt, nF_copy);
       nF_copy++;
    }
    
    int faceID = PUMI_MeshEnt_ID(meshEnt);
// get vertices from adjacency
    std::vector<pMeshEnt> vecVtx;
    PUMI_MeshEnt_GetAdj(meshEnt, PUMI_VERTEX, 0, vecVtx);
    int iNumVtx = vecVtx.size();
    for (int iVtx = 0; iVtx < iNumVtx; ++iVtx) {
      int VtxID = PUMI_MeshEnt_ID(vecVtx[iVtx]);
      mesh.elementBoundaryNodesArray[faceID*mesh.nNodes_elementBoundary+iVtx] = VtxID;  //populating the element boundary (face) nodes array (3 is hardcode for tets)
    }

//get regions from adjacency
    std::vector<pMeshEnt> vecRgn;
    PUMI_MeshEnt_GetAdj(meshEnt, PUMI_REGION, 0, vecRgn);
    int iNumRgn = vecRgn.size();
    int RgnID[2] = {-1,-1}; 
    int LocalFaceNumber[2] = {-1,-1};
    for (int iRgn = 0; iRgn < iNumRgn; ++iRgn) {
      RgnID[iRgn] = PUMI_MeshEnt_ID(vecRgn[iRgn]);
      mesh.elementBoundaryElementsArray[faceID*2+iRgn]= RgnID[iRgn];

      //following is tricky. it stores local face number of the face we are already on for this element (vecRgn[iVtx])
      std::vector<pMeshEnt> vecFace;
      PUMI_MeshEnt_GetAdj(vecRgn[iRgn], PUMI_FACE, 0, vecFace);
      int iNumFace = vecFace.size();
      for (int iFace = 0; iFace < iNumFace; ++iFace) {
         if(vecFace[iFace] == meshEnt) {//if face found
//             LocalFaceNumber[iRgn] = iFace; //old method caused conflict in face numbers
//following is going to look adhoc and arbitrary but it is needed to resolve a conflict between the way proteus handles faces and we handle faces
//The order in which we retrieve faces from adjacency is different to what proteus expects them to be in their elementBoundaries array. To resolve this
//we need to change the local number of the face of an element to what proteus expects it to be and then it works. Whussh. magic.             
             if(iFace==0) LocalFaceNumber[iRgn]=3;
             if(iFace==1) LocalFaceNumber[iRgn]=2;
             if(iFace==2) LocalFaceNumber[iRgn]=0;
             if(iFace==3) LocalFaceNumber[iRgn]=1;
         }
      }
      assert(LocalFaceNumber[iRgn] != -1);
      mesh.elementBoundaryLocalElementBoundariesArray[faceID*2+iRgn]=LocalFaceNumber[iRgn]; 
    }

    //left and right regions are shared by this face we are currntly on
    int leftRgnID = RgnID[0]; int leftLocalFaceNumber = LocalFaceNumber[0];
    int rightRgnID = RgnID[1]; int rightLocalFaceNumber = LocalFaceNumber[1];

    //left region is always there, so either rightRgnID will have an actual ID if this face is shared, or will contain -1 if it is an exterior face
    mesh.elementNeighborsArray[leftRgnID*mesh.nElementBoundaries_element+leftLocalFaceNumber] = rightRgnID; 
    mesh.elementBoundariesArray[leftRgnID*mesh.nElementBoundaries_element+leftLocalFaceNumber] = faceID;

    //if only 1 region is adjacent to this face, that means it is an exterior face (or on part boundary, but parallel is for later) todo: parallel
    if(iNumRgn==1 && RgnID[1]==-1 && LocalFaceNumber[1]==-1) { //last 2 checks are only for sanity
      mesh.elementBoundaryElementsArray[faceID*2+1] = -1; 
      mesh.elementBoundaryLocalElementBoundariesArray[faceID*2+1] = -1;
      exteriorElementBoundaries.insert(faceID); //exterior face as only 1 region adjacent
    } else { //2 regions are shared by this face so interior face
      mesh.elementNeighborsArray[rightRgnID*mesh.nElementBoundaries_element+rightLocalFaceNumber] = leftRgnID;
      mesh.elementBoundariesArray[rightRgnID*mesh.nElementBoundaries_element+rightLocalFaceNumber] = faceID;
      interiorElementBoundaries.insert(faceID);
    }

    PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
    nF++;
  }
  PUMI_PartEntIter_Del(EntIt);
  assert(nF_owned=faces_owned-1);
  assert(nF_copy=mesh.nElementBoundaries_global-1);
 
  //construct interior and exterior element boundaries array 
  mesh.nInteriorElementBoundaries_global = interiorElementBoundaries.size();
  mesh.interiorElementBoundariesArray = new int[mesh.nInteriorElementBoundaries_global];
  mesh.nExteriorElementBoundaries_global = exteriorElementBoundaries.size();
  mesh.exteriorElementBoundariesArray = new int[mesh.nExteriorElementBoundaries_global];

  int ebNI=0,ebNE=0;
  for (std::set<int>::iterator ebN=interiorElementBoundaries.begin();ebN != interiorElementBoundaries.end(); ebN++,ebNI++)
      mesh.interiorElementBoundariesArray[ebNI] = *ebN;
  for (std::set<int>::iterator ebN=exteriorElementBoundaries.begin();ebN != exteriorElementBoundaries.end(); ebN++,ebNE++)
      mesh.exteriorElementBoundariesArray[ebNE] = *ebN;
  
  return 0;
//////end: build face list (elementBoundary and nodeBoundary arrays) 
//
}

int MeshAdaptPUMIDrvr::ConstructEdges(Mesh& mesh) {
/////////build edge data structutes and interior and exterior element boundary array

  std::cout << "Constructing edge data structures\n" ; 
  mesh.edgeNodesArray = new int[mesh.nEdges_global*2];

  pPartEntIter EntIt;
  pMeshEnt meshEnt;
  PUMI_PartEntIter_Init (PUMI_Part, PUMI_EDGE, PUMI_ALLTOPO, EntIt);
  int isEnd = 0;
  int nE = 0; int nE_owned = 0;
  int nE_copy = 0;
  while (!isEnd) {
     PUMI_PartEntIter_GetNext(EntIt, meshEnt);
     
     int owned;
     PUMI_MeshEnt_IsOwned(meshEnt, PUMI_Part, &owned);
     if(owned) {
       PUMI_MeshEnt_SetID(meshEnt, nE_owned);
       nE_owned++;
     } else {
       PUMI_MeshEnt_SetID(meshEnt, nE_copy);
       nE_copy++;
     }
     
     int edgeID = PUMI_MeshEnt_ID(meshEnt);
     std::vector<pMeshEnt> vecVtx;
     PUMI_MeshEnt_GetAdj(meshEnt, PUMI_VERTEX, 0, vecVtx);
     int iNumVtx = vecVtx.size();
     for(int iVtx=0; iVtx<iNumVtx; iVtx++) {
        int VtxID = PUMI_MeshEnt_ID(vecVtx[iVtx]);
        mesh.edgeNodesArray[edgeID*2+iVtx] = VtxID; 
     }

     PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
     nE++;
  }
  PUMI_PartEntIter_Del(EntIt);
  
//  std::cout << "\n" << nE << "\n";
  std::vector<std::set<int> > nodeStar(mesh.nNodes_global);
  for (int edgeN=0;edgeN<mesh.nEdges_global;edgeN++) {
//      std::cout << edgeN << " " << mesh.edgeNodesArray[edgeN*2+0] << " " << mesh.edgeNodesArray[edgeN*2+1] << "\n";
      nodeStar[mesh.edgeNodesArray[edgeN*2+0]].insert(mesh.edgeNodesArray[edgeN*2+1]);
      nodeStar[mesh.edgeNodesArray[edgeN*2+1]].insert(mesh.edgeNodesArray[edgeN*2+0]);
  }

//not completely comfortable with folllowing implementation
  mesh.nodeStarOffsets = new int[mesh.nNodes_global+1];
  mesh.nodeStarOffsets[0] = 0;
  for (int nN=1;nN<mesh.nNodes_global+1;nN++) {
     mesh.nodeStarOffsets[nN] = mesh.nodeStarOffsets[nN-1] + nodeStar[nN-1].size();
//     std::cout << "nodeStarOffesets: " << nN << ": " << mesh.nodeStarOffsets[nN] << "\n"; 
  }
  
  mesh.nodeStarArray = new int[mesh.nodeStarOffsets[mesh.nNodes_global]];
  for (int nN=0,offset=0;nN<mesh.nNodes_global;nN++) {
      for (std::set<int>::iterator nN_star=nodeStar[nN].begin();nN_star!=nodeStar[nN].end();nN_star++,offset++) {
          mesh.nodeStarArray[offset] = *nN_star;
//          std::cout << "nodeStarArray: " << offset << ": " << mesh.nodeStarArray[offset] << std::endl;
      }
  }

  mesh.max_nNodeNeighbors_node=0;
  for (int nN=0;nN<mesh.nNodes_global;nN++)
      mesh.max_nNodeNeighbors_node=std::max(mesh.max_nNodeNeighbors_node,mesh.nodeStarOffsets[nN+1]-mesh.nodeStarOffsets[nN]);
  //mwf repeat for node-->elements arrays
  std::vector<std::set<int> > nodeElementsStar(mesh.nNodes_global);
  for (int eN = 0; eN < mesh.nElements_global; eN++)  {
     for (int nN = 0; nN < mesh.nNodes_element; nN++) {
//         std::cerr << eN << " " << nN << " " << mesh.elementNodesArray[eN*mesh.nNodes_element+nN] << "\n";
         nodeElementsStar[mesh.elementNodesArray[eN*mesh.nNodes_element+nN]].insert(eN);
     }
  }
  mesh.nodeElementOffsets = new int[mesh.nNodes_global+1];
  mesh.nodeElementOffsets[0] = 0;
  for (int nN = 0; nN < mesh.nNodes_global; nN++)
     mesh.nodeElementOffsets[nN+1] = mesh.nodeElementOffsets[nN]+nodeElementsStar[nN].size();
  mesh.nodeElementsArray  = new int[mesh.nodeElementOffsets[mesh.nNodes_global]];
  for (int nN=0,offset=0; nN < mesh.nNodes_global; nN++)   {
     for (std::set<int>::iterator eN_star = nodeElementsStar[nN].begin(); eN_star != nodeElementsStar[nN].end(); eN_star++,offset++)   {
          mesh.nodeElementsArray[offset] = *eN_star;
     }
  }
  //mwf end node-->elements construction
  
  assert(nE_owned=edges_owned-1);
  assert(nE_copy=mesh.nEdges_global-1);
  return 0;
//end: build edge data structures (also added star data structures here itself)
}

int MeshAdaptPUMIDrvr::ConstructMaterialArrays(Mesh& mesh) {

  std::cout << "Constructing material data structures\n" ; 

  mesh.elementBoundaryMaterialTypes = new int[mesh.nElementBoundaries_global];
  mesh.nodeMaterialTypes = new int[mesh.nNodes_global];
  mesh.elementMaterialTypes = new int[mesh.nElements_global];

  memset(mesh.nodeMaterialTypes,DEFAULT_NODE_MATERIAL,mesh.nNodes_global*sizeof(int));
  memset(mesh.elementMaterialTypes,DEFAULT_ELEMENT_MATERIAL,mesh.nElements_global*sizeof(int));
 
  pPartEntIter EntIt;
  pMeshEnt meshEnt;
  PUMI_PartEntIter_Init (PUMI_Part, PUMI_FACE, PUMI_ALLTOPO, EntIt);
  int isEnd = 0;
  int nF = 0;
  //populate elementBoundary Material arrays
  while (!isEnd) {
     PUMI_PartEntIter_GetNext(EntIt, meshEnt);
     
     int faceID = PUMI_MeshEnt_ID(meshEnt); 
     int geomType; 
     PUMI_MeshEnt_GetGeomClasType (meshEnt, &geomType);
        
     if(geomType==2)
        mesh.elementBoundaryMaterialTypes[faceID] = EXTERIOR_ELEMENT_BOUNDARY_MATERIAL;
     if(geomType==3)
        mesh.elementBoundaryMaterialTypes[faceID] = INTERIOR_ELEMENT_BOUNDARY_MATERIAL;
     PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
     nF++;
  }
  PUMI_PartEntIter_Del(EntIt);
  
  //populate node Material arrays
  PUMI_PartEntIter_Init (PUMI_Part, PUMI_VERTEX, PUMI_ALLTOPO, EntIt);
  isEnd = 0;
  int nN = 0;
  //populate elementBoundary Material arrays
  while (!isEnd) {
     PUMI_PartEntIter_GetNext(EntIt, meshEnt);

     int vtxID = PUMI_MeshEnt_ID(meshEnt);
     int geomType;
     PUMI_MeshEnt_GetGeomClasType (meshEnt, &geomType);

     if(geomType==2 || geomType==1 || geomType==0)
        mesh.nodeMaterialTypes[vtxID] = EXTERIOR_NODE_MATERIAL;
     if(geomType==3)
        mesh.nodeMaterialTypes[vtxID] = INTERIOR_NODE_MATERIAL;
     PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
     nN++;
  }
  PUMI_PartEntIter_Del(EntIt);

  return 0;

}




