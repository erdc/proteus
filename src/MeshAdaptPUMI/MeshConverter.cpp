#include <algorithm>

#include "MeshAdaptPUMI.h"
#include "mesh.h"

const int DEFAULT_ELEMENT_MATERIAL=0;
const int DEFAULT_NODE_MATERIAL=-1;
const int INTERIOR_NODE_MATERIAL=0;
const int EXTERIOR_NODE_MATERIAL=1;
const int INTERIOR_ELEMENT_BOUNDARY_MATERIAL=0;
const int EXTERIOR_ELEMENT_BOUNDARY_MATERIAL=1;

//Main API to construct a serial pumi mesh, what we do is contruct the global mesh when working with serial
//and let Proteus populate the subdomain data structures (though it will be exactly same)
int MeshAdaptPUMIDrvr::ConstructFromSerialPUMIMesh(Mesh& mesh)
{

  std::cout << "Constructing global data structures\n"; 
  int numGlobElem = m->count(3);
  mesh.nElements_global = numGlobElem;
  elms_owned = numGlobElem;

  int numGlobNodes = m->count(0);
  mesh.nNodes_global = numGlobNodes;
  vtx_owned = numGlobNodes;

  int numGlobFaces = m->count(2);
  mesh.nElementBoundaries_global = numGlobFaces;
  faces_owned = numGlobFaces;

  int numGlobEdges = m->count(1);
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

//Construct node data structures for Proteus
int MeshAdaptPUMIDrvr::ConstructNodes(Mesh& mesh) {
//////build node array (coordinates?)

//set interior and exterior node material flags after get boundary info in
//constructElementBoundaryElementsArray_edge 
//if nodeMaterialTypes left as DEFAULT_NODE_MATERIALi
  mesh.nodeArray = new double[mesh.nNodes_global*3];
  memset(mesh.nodeArray,0,mesh.nNodes_global*3*sizeof(double));
  mesh.nodeMaterialTypes = new int[mesh.nNodes_global];
  memset(mesh.nodeMaterialTypes,DEFAULT_NODE_MATERIAL,mesh.nNodes_global*sizeof(int));

  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;
  size_t i = 0;
  while ((e = m->iterate(it))) {
    apf::Vector3 x;
    m->getPoint(e, 0, x);
    for(int j=0; j<3; j++)
      mesh.nodeArray[i * 3 + j]= x[j];
    ++i;
  }
  m->end(it);

  freeNumbering(local[0]);
  local[0] = apf::numberOverlapNodes(m, "proteus_local_0");
  
  return 0;
} 

int MeshAdaptPUMIDrvr::ConstructElements(Mesh& mesh)
{
//////build node list for the elements
  mesh.elementNodesArray = new int[mesh.nElements_global*mesh.nNodes_element];
  mesh.elementMaterialTypes = new int[mesh.nElements_global];
  memset(mesh.elementMaterialTypes,DEFAULT_ELEMENT_MATERIAL,mesh.nElements_global*sizeof(int));

  apf::MeshIterator* it = m->begin(m->getDimension());
  apf::MeshEntity* e;
  size_t i = 0;
  while ((e = m->iterate(it))) {
    apf::Downward v;
    int iNumVtx = m->getDownward(e, 0, v);
    for (int j = 0; j < iNumVtx; ++j) {
      int vtxID = apf::getNumber(local[0], v[j], 0, 0);
      mesh.elementNodesArray[i * mesh.nNodes_element + j] = vtxID;
    }
    ++i;
  } //region loop
  m->end(it);

  freeNumbering(local[3]);
  local[3] = apf::numberElements(m, "proteus_local_3");

  return 0;
}

/* following is going to look adhoc and arbitrary but it is needed to
   resolve a conflict between the way proteus handles faces and we handle faces.
   The order in which we retrieve faces from adjacency is different to what
   proteus expects them to be in their elementBoundaries array.
   To resolve this we need to change the local number of the face of an element
   to what proteus expects it to be and then it works. Whussh. magic.

   original comment above preserved for entertainment.
   This maps SCOREC's tet face numbering to that of proteus.
 */

int getProteusFaceIdx(apf::Mesh* m, apf::MeshEntity* e, apf::MeshEntity* f)
{
  apf::Downward fs;
  int nfs = m->getDownward(e, 2, fs);
  int idx_apf = apf::findIn(fs, nfs, f);
  static int const map[4] = {3,2,0,1};
  return map[idx_apf];
}

int MeshAdaptPUMIDrvr::ConstructBoundaries(Mesh& mesh)
{
//build face list (elementBoundary and nodeBoundary arrays)
//this is a bit complex and should be verified
//Enter at your own peril for those who stray will be lost
  std::set<int> interiorElementBoundaries;
 //do we need this? we can do this by looping over model faces, right?
  std::set<int> exteriorElementBoundaries;

  mesh.elementBoundaryNodesArray =
    new int[mesh.nElementBoundaries_global * mesh.nNodes_elementBoundary];
  mesh.elementBoundaryElementsArray =
    new int[mesh.nElementBoundaries_global *2];
  mesh.elementBoundaryLocalElementBoundariesArray =
    new int[mesh.nElementBoundaries_global*2];
  mesh.elementNeighborsArray =
    new int[mesh.nElements_global * mesh.nElementBoundaries_element];
  mesh.elementBoundariesArray =
    new int[mesh.nElements_global * mesh.nElementBoundaries_element];

  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* f;
  size_t i = 0;
  while ((f = m->iterate(it))) {
// get vertices from adjacency
    apf::Downward vs;
    int iNumVtx = m->getDownward(f, 0, vs);
    for (int iVtx = 0; iVtx < iNumVtx; ++iVtx) {
      int vtxID = apf::getNumber(local[0], vs[iVtx], 0, 0);
      mesh.elementBoundaryNodesArray[
        i * mesh.nNodes_elementBoundary + iVtx] = vtxID;
    }
//get regions from adjacency
    apf::Up rs;
    m->getUp(f, rs);
    int iNumRgn = rs.n;
    int RgnID[2] = {-1,-1}; 
    int LocalFaceNumber[2] = {-1,-1};
    for (int iRgn = 0; iRgn < iNumRgn; ++iRgn) {
      RgnID[iRgn] = apf::getNumber(local[3], rs.e[iRgn], 0, 0);
      mesh.elementBoundaryElementsArray[i * 2 + iRgn]= RgnID[iRgn];
      LocalFaceNumber[iRgn] = getProteusFaceIdx(m, rs.e[iRgn], f);
      assert(LocalFaceNumber[iRgn]!=-1);
      mesh.elementBoundaryLocalElementBoundariesArray[
        i * 2 + iRgn] = LocalFaceNumber[iRgn];
    }
    //left and right regions are shared by this face we are currntly on
    int leftRgnID = RgnID[0]; int leftLocalFaceNumber = LocalFaceNumber[0];
    int rightRgnID = RgnID[1]; int rightLocalFaceNumber = LocalFaceNumber[1];

    /*left region is always there, so either rightRgnID will have
      an actual ID if this face is shared, or will contain -1
      if it is an exterior face */
    mesh.elementNeighborsArray[
      leftRgnID * mesh.nElementBoundaries_element + leftLocalFaceNumber]
      = rightRgnID;
    mesh.elementBoundariesArray[
      leftRgnID * mesh.nElementBoundaries_element + leftLocalFaceNumber]
      = i;

    /* if only 1 region is adjacent to this face,
       that means it is an exterior face */
    if(iNumRgn==1) {
      assert(RgnID[1]==-1);
      assert(LocalFaceNumber[1]==-1); //last 2 checks are only for sanity
      mesh.elementBoundaryElementsArray[i * 2 + 1] = -1;
      mesh.elementBoundaryLocalElementBoundariesArray[i * 2 + 1] = -1;
      //exterior face as only 1 region adjacent
      exteriorElementBoundaries.insert(i);
    } else { //2 regions are shared by this face so interior face
      mesh.elementNeighborsArray[
        rightRgnID * mesh.nElementBoundaries_element + rightLocalFaceNumber]
        = leftRgnID;
      mesh.elementBoundariesArray[
        rightRgnID * mesh.nElementBoundaries_element + rightLocalFaceNumber]
        = i;
      interiorElementBoundaries.insert(i);
    }
    i++;
  }
  m->end(it);
 
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
}

/* these algorithms are totally independent of SCOREC;
   they form the proteus
   vertex to vertex and vertex to element adjacency tables
   from the edge to vertex and element to vertex tables */
static void createStars(Mesh& mesh)
{
  std::vector<std::set<int> > nodeStar(mesh.nNodes_global);
  for (int edgeN = 0; edgeN < mesh.nEdges_global; edgeN++) {
    nodeStar[mesh.edgeNodesArray[edgeN * 2 + 0]].insert(
             mesh.edgeNodesArray[edgeN * 2 + 1]);
    nodeStar[mesh.edgeNodesArray[edgeN * 2 + 1]].insert(
             mesh.edgeNodesArray[edgeN * 2 + 0]);
  }

  mesh.nodeStarOffsets = new int[mesh.nNodes_global + 1];
  mesh.nodeStarOffsets[0] = 0;
  for (int nN = 1; nN <= mesh.nNodes_global; nN++)
    mesh.nodeStarOffsets[nN] =
      mesh.nodeStarOffsets[nN - 1] + nodeStar[nN - 1].size();
  
  mesh.nodeStarArray = new int[mesh.nodeStarOffsets[mesh.nNodes_global]];
  for (int nN = 0, offset = 0; nN < mesh.nNodes_global; nN++) {
    for (std::set<int>::iterator nN_star=nodeStar[nN].begin();
         nN_star!=nodeStar[nN].end();
         nN_star++, offset++)
       mesh.nodeStarArray[offset] = *nN_star;
  }

  mesh.max_nNodeNeighbors_node = 0;
  for (int nN = 0; nN < mesh.nNodes_global; nN++)
      mesh.max_nNodeNeighbors_node =
        std::max(mesh.max_nNodeNeighbors_node,
                 mesh.nodeStarOffsets[nN + 1] - mesh.nodeStarOffsets[nN]);

  std::vector<std::set<int> > nodeElementsStar(mesh.nNodes_global);
  for (int eN = 0; eN < mesh.nElements_global; eN++) {
    for (int nN = 0; nN < mesh.nNodes_element; nN++) {
      nodeElementsStar[
        mesh.elementNodesArray[eN * mesh.nNodes_element + nN]
        ].insert(eN);
    }
  }
  mesh.nodeElementOffsets = new int[mesh.nNodes_global + 1];
  mesh.nodeElementOffsets[0] = 0;
  for (int nN = 0; nN < mesh.nNodes_global; nN++)
     mesh.nodeElementOffsets[nN + 1] =
       mesh.nodeElementOffsets[nN] + nodeElementsStar[nN].size();
  mesh.nodeElementsArray = new int[mesh.nodeElementOffsets[mesh.nNodes_global]];
  for (int nN=0,offset=0; nN < mesh.nNodes_global; nN++)   {
    for (std::set<int>::iterator eN_star = nodeElementsStar[nN].begin();
         eN_star != nodeElementsStar[nN].end();
         eN_star++, offset++)
      mesh.nodeElementsArray[offset] = *eN_star;
  }
}

int MeshAdaptPUMIDrvr::ConstructEdges(Mesh& mesh)
{
/////////build edge data structutes and interior and exterior element boundary array
  mesh.edgeNodesArray = new int[mesh.nEdges_global * 2];

  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* e;
  size_t i = 0;
  while ((e = m->iterate(it))) {
    apf::MeshEntity* v[2];
    m->getDownward(e, 0, v);
    for(int iVtx=0; iVtx < 2; ++iVtx) {
      int vtxID = apf::getNumber(local[0], v[iVtx], 0, 0);
      mesh.edgeNodesArray[i * 2 + iVtx] = vtxID;
    }
    i++;
  }
  m->end(it);
  
//end: build edge data structures (also added star data structures here itself)
  createStars(mesh);
  
  return 0;
}

int MeshAdaptPUMIDrvr::ConstructMaterialArrays(Mesh& mesh)
{
  mesh.elementBoundaryMaterialTypes = new int[mesh.nElementBoundaries_global];
  mesh.nodeMaterialTypes = new int[mesh.nNodes_global];
  mesh.elementMaterialTypes = new int[mesh.nElements_global];

  memset(mesh.nodeMaterialTypes,DEFAULT_NODE_MATERIAL,
      mesh.nNodes_global*sizeof(int));
  memset(mesh.elementMaterialTypes,DEFAULT_ELEMENT_MATERIAL,
      mesh.nElements_global*sizeof(int));
 
  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* f;
  //populate elementBoundary Material arrays
  static int const material_table[4] = {
    EXTERIOR_ELEMENT_BOUNDARY_MATERIAL,
    EXTERIOR_ELEMENT_BOUNDARY_MATERIAL,
    EXTERIOR_ELEMENT_BOUNDARY_MATERIAL,
    INTERIOR_ELEMENT_BOUNDARY_MATERIAL
  };
  size_t i = 0;
  while ((f = m->iterate(it))) {
    int geomType = m->getModelType(m->toModel(f));
    mesh.elementBoundaryMaterialTypes[i] = material_table[geomType];
    i++;
  }
  m->end(it);
  
  //populate node Material arrays
  it = m->begin(0);
  apf::MeshEntity* v;
  //populate elementBoundary Material arrays
  i = 0;
  while ((v = m->iterate(it))) {
    int geomType = m->getModelType(m->toModel(v));
    mesh.nodeMaterialTypes[i] = material_table[geomType];
    i++;
  }
  m->end(it);

  return 0;
}

//function to update the material arrays for boundary conditions
int MeshAdaptPUMIDrvr::UpdateMaterialArrays(Mesh& mesh, int bdryID, int geomTag)
{
  apf::ModelEntity* geomEnt = m->findModelEntity(2, geomTag);
  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* f;
  //populate elementBoundary Material arrays
  size_t i = 0;
  while ((f = m->iterate(it))) {
    if (m->toModel(f) == geomEnt) {
      mesh.elementBoundaryMaterialTypes[i] = bdryID;
      apf::Downward vs;
      int iNumVtx = m->getDownward(f, 0, vs);
      for(int iVtx=0; iVtx < iNumVtx; ++iVtx) {
        int vtxID = apf::getNumber(local[0], vs[iVtx], 0, 0);
        mesh.nodeMaterialTypes[vtxID] = bdryID;
      }
    }
    ++i;
  }
  m->end(it);

  return 0;
}
