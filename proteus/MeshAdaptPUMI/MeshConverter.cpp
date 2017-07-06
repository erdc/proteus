#include <algorithm>

#include "MeshAdaptPUMI.h"
#include "mesh.h"
#include <apfShape.h>

#include <sstream>

static apf::Numbering* numberOwnedEntitiesFirst(apf::Mesh* m, int dimension,int initialReconstructed)
{
  std::stringstream ss;
  ss << "proteus_number_" << dimension;
  std::string s = ss.str();
  apf::FieldShape* shape;
  if (dimension) /* this switch is just to help rendering */
    shape = apf::getConstant(dimension);
  else
    shape = m->getShape();
  apf::Numbering* n = createNumbering(m, s.c_str(), shape, 1);
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(dimension);
  int i = 0;
  if(initialReconstructed){
    while ((e = m->iterate(it)))
      apf::number(n, e, 0, 0, i++);
    m->end(it);
  }
  else{
    while ((e = m->iterate(it)))
      if (m->isOwned(e))
        apf::number(n, e, 0, 0, i++);
    m->end(it);
    it = m->begin(dimension);
    while ((e = m->iterate(it)))
      if (!m->isOwned(e))
        apf::number(n, e, 0, 0, i++);
    m->end(it);
  }
  return n;
}

//Main API to construct a serial pumi mesh,
//what we do is contruct the global mesh when working with serial
//and let Proteus populate the subdomain data structures
//(though it will be exactly same)
int MeshAdaptPUMIDrvr::constructFromSerialPUMIMesh(Mesh& mesh)
{
  assert(m != 0);
  std::cout << "Constructing global data structures\n"; 

  int dim = m->getDimension();
  mesh.nElements_global = m->count(dim);

  mesh.nNodes_global = m->count(0);

  mesh.nElementBoundaries_global = m->count(dim - 1);

  mesh.nEdges_global = m->count(1);

//nNodes_element for now is constant for the entire mesh, Ask proteus about using mixed meshes
  switch (dim) {
    case 2:
      mesh.nNodes_element = 3;
      mesh.nNodes_elementBoundary = 2;
      mesh.nElementBoundaries_element = 3;
      break;
    case 3:
      mesh.nNodes_element = 4;
      mesh.nNodes_elementBoundary = 3;
      mesh.nElementBoundaries_element = 4;
      break;
    default:
      apf::fail("dimension is not 2 or 3\n");
      break;
  }
#ifdef MESH_INFO
  std::cerr << "*******Proteus Mesh Stats*********\n";
  std::cerr << "Number of elements " << mesh.nElements_global << "\n";
  std::cerr << "Number of nodes " << mesh.nNodes_global << "\n";
  std::cerr << "Number of boundaries " << mesh.nElementBoundaries_global << "\n";
  std::cerr << "Number of edges " << mesh.nEdges_global << "\n";
#endif
  
  numberLocally();
  constructNodes(mesh);
  constructElements(mesh);
  constructBoundaries(mesh);
  constructEdges(mesh);
  constructMaterialArrays(mesh);

  return 0;
}

void MeshAdaptPUMIDrvr::numberLocally()
{
  for (int d = 0; d <= m->getDimension(); ++d) {
    freeNumbering(local[d]);
    local[d] = numberOwnedEntitiesFirst(m, d,initialReconstructed);
  }
  if(initialReconstructed) 
    initialReconstructed = 0;
}

int MeshAdaptPUMIDrvr::localNumber(apf::MeshEntity* e)
{
  return getNumber(local[apf::getDimension(m, e)], e, 0, 0);
}

int MeshAdaptPUMIDrvr::constructNodes(Mesh& mesh)
{
  mesh.nodeArray = new double[mesh.nNodes_global * 3];
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    int i = localNumber(e);
    apf::Vector3 x;
    m->getPoint(e, 0, x);
    for(int j=0; j<3; j++)
      mesh.nodeArray[i * 3 + j]= x[j];
  }
  m->end(it);
  return 0;
} 

int MeshAdaptPUMIDrvr::constructElements(Mesh& mesh)
{
  mesh.elementNodesArray = new int[mesh.nElements_global*mesh.nNodes_element];
  apf::MeshIterator* it = m->begin(m->getDimension());
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    int i = localNumber(e);
    apf::Downward v;
    int iNumVtx = m->getDownward(e, 0, v);
    for (int j = 0; j < iNumVtx; ++j) {
      int vtxID = localNumber(v[j]);
      mesh.elementNodesArray[i * mesh.nNodes_element + j] = vtxID;
    }
  }
  m->end(it);
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

int getProteusBoundaryIdx(apf::Mesh* m, apf::MeshEntity* e, apf::MeshEntity* f)
{
  apf::Downward fs;
  int dim = m->getDimension();
  int nfs = m->getDownward(e, dim - 1, fs);
  int idx_apf = apf::findIn(fs, nfs, f);
  /* Proteus convention is that the face index equals the vertex index
     of the vertex opposite to the face.
     Proteus and PUMI should have consistent vertex orderings for
     simplices, but the above rule makes the side orderings different */
  static int const tet_boundary_map[4] = {3,2,0,1};
  static int const tri_boundary_map[3] = {2,0,1};
  static int const* const boundary_maps[4] = {
    0,
    0,
    tri_boundary_map,
    tet_boundary_map
  };
  return boundary_maps[dim][idx_apf];
}

int MeshAdaptPUMIDrvr::constructBoundaries(Mesh& mesh)
{
//build face list (elementBoundary and nodeBoundary arrays)
//Enter at your own peril for those who stray will be lost
  std::set<int> interiorElementBoundaries;
  std::set<int> exteriorElementBoundaries;

  mesh.elementBoundaryNodesArray =
    new int[mesh.nElementBoundaries_global * mesh.nNodes_elementBoundary];
  mesh.elementBoundaryElementsArray =
    new int[mesh.nElementBoundaries_global * 2];
  mesh.elementBoundaryLocalElementBoundariesArray =
    new int[mesh.nElementBoundaries_global * 2];
  mesh.elementNeighborsArray =
    new int[mesh.nElements_global * mesh.nElementBoundaries_element];
  mesh.elementBoundariesArray =
    new int[mesh.nElements_global * mesh.nElementBoundaries_element];
  exteriorGlobaltoLocalElementBoundariesArray = 
    new int[mesh.nElementBoundaries_global];
  int exterior_count = 0; //counter for external boundaries

  int dim = m->getDimension();
  apf::MeshIterator* it = m->begin(dim - 1);
  apf::MeshEntity* f;
  while ((f = m->iterate(it))) {
    int i = localNumber(f);
// get vertices from adjacency
    apf::Downward vs;
    int iNumVtx = m->getDownward(f, 0, vs);
    for (int iVtx = 0; iVtx < iNumVtx; ++iVtx) {
      int vtxID = localNumber(vs[iVtx]);
      mesh.elementBoundaryNodesArray[
        i * mesh.nNodes_elementBoundary + iVtx] = vtxID;
    }
//get regions from adjacency
    apf::Up rs;
    m->getUp(f, rs);
    int iNumRgn = rs.n;
    int RgnID[2] = {-1,-1}; 
    int localBoundaryNumber[2] = {-1,-1};
    for (int iRgn = 0; iRgn < iNumRgn; ++iRgn) {
      RgnID[iRgn] = localNumber(rs.e[iRgn]);
      mesh.elementBoundaryElementsArray[i * 2 + iRgn]= RgnID[iRgn];
      localBoundaryNumber[iRgn] = getProteusBoundaryIdx(m, rs.e[iRgn], f);
      assert(localBoundaryNumber[iRgn] != -1);
      mesh.elementBoundaryLocalElementBoundariesArray[
        i * 2 + iRgn] = localBoundaryNumber[iRgn];
    }
    //left and right regions are shared by this face we are currntly on
    int leftRgnID = RgnID[0]; int leftLocalBoundaryNumber = localBoundaryNumber[0];
    int rightRgnID = RgnID[1]; int rightLocalBoundaryNumber = localBoundaryNumber[1];
    /*left region is always there, so either rightRgnID will have
      an actual ID if this face is shared, or will contain -1
      if it is an exterior face */
    mesh.elementNeighborsArray[
      leftRgnID * mesh.nElementBoundaries_element + leftLocalBoundaryNumber]
      = rightRgnID;
    mesh.elementBoundariesArray[
      leftRgnID * mesh.nElementBoundaries_element + leftLocalBoundaryNumber]
      = i;

    /* if only 1 region is adjacent to this face,
       that means it is an exterior face */
    if(iNumRgn==1) {
      assert(RgnID[1]==-1);
      assert(localBoundaryNumber[1]==-1); //last 2 checks are only for sanity
      mesh.elementBoundaryElementsArray[i * 2 + 1] = -1;
      mesh.elementBoundaryLocalElementBoundariesArray[i * 2 + 1] = -1;
      //exterior face as only 1 region adjacent
      exteriorElementBoundaries.insert(i);

      //construct inverse mapping 
      exteriorGlobaltoLocalElementBoundariesArray[i] = exterior_count;
      exterior_count++;

    } else { //2 regions are shared by this face so interior face
      mesh.elementNeighborsArray[
        rightRgnID * mesh.nElementBoundaries_element + rightLocalBoundaryNumber]
        = leftRgnID;
      mesh.elementBoundariesArray[
        rightRgnID * mesh.nElementBoundaries_element + rightLocalBoundaryNumber]
        = i;
      interiorElementBoundaries.insert(i);
    }
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
  for (int nN = 0, offset = 0; nN < mesh.nNodes_global; nN++)
    for (std::set<int>::iterator nN_star=nodeStar[nN].begin();
         nN_star!=nodeStar[nN].end();
         nN_star++, offset++)
       mesh.nodeStarArray[offset] = *nN_star;

  mesh.max_nNodeNeighbors_node = 0;
  for (int nN = 0; nN < mesh.nNodes_global; nN++)
      mesh.max_nNodeNeighbors_node =
          std::max(mesh.max_nNodeNeighbors_node,
                   mesh.nodeStarOffsets[nN + 1] - mesh.nodeStarOffsets[nN]);

  std::vector<std::set<int> > nodeElementsStar(mesh.nNodes_global);
  for (int eN = 0; eN < mesh.nElements_global; eN++)
    for (int nN = 0; nN < mesh.nNodes_element; nN++)
      nodeElementsStar[
          mesh.elementNodesArray[eN * mesh.nNodes_element + nN]
          ].insert(eN);
  mesh.nodeElementOffsets = new int[mesh.nNodes_global + 1];
  mesh.nodeElementOffsets[0] = 0;
  for (int nN = 0; nN < mesh.nNodes_global; nN++)
     mesh.nodeElementOffsets[nN + 1] =
       mesh.nodeElementOffsets[nN] + nodeElementsStar[nN].size();
  mesh.nodeElementsArray = new int[mesh.nodeElementOffsets[mesh.nNodes_global]];
  for (int nN=0,offset=0; nN < mesh.nNodes_global; nN++)
    for (std::set<int>::iterator eN_star = nodeElementsStar[nN].begin();
         eN_star != nodeElementsStar[nN].end();
         eN_star++, offset++)
      mesh.nodeElementsArray[offset] = *eN_star;
}

int MeshAdaptPUMIDrvr::constructEdges(Mesh& mesh)
{
  mesh.edgeNodesArray = new int[mesh.nEdges_global * 2];
  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    int i = localNumber(e);
    apf::MeshEntity* v[2];
    m->getDownward(e, 0, v);
    for(int iVtx=0; iVtx < 2; ++iVtx) {
      int vtxID = localNumber(v[iVtx]);
      mesh.edgeNodesArray[i * 2 + iVtx] = vtxID;
    }
  }
  m->end(it);
  createStars(mesh);
  return 0;
}

#define INTERIOR_MATERIAL 0
#define EXTERIOR_MATERIAL 1
#define DEFAULT_ELEMENT_MATERIAL INTERIOR_MATERIAL

static int getInOutMaterial(apf::Mesh* m, apf::MeshEntity* e)
{
  if (m->getModelType(m->toModel(e)) == m->getDimension())
    return INTERIOR_MATERIAL;
  else
    return EXTERIOR_MATERIAL;
}

/* This function builds the element, elementBoundary, and node
 * Material arrays and fills them in with zero or one depending
 * one whether the entity is classified on the geometric boundary
 * or not.
 * This forms a baseline material tagging for all entities,
 * which later gets overwritten for some entities with more
 * specific values in MeshAdaptPUMIDrvr::UpdateMaterialArrays.
 */
int MeshAdaptPUMIDrvr::constructMaterialArrays(Mesh& mesh)
{
  mesh.elementBoundaryMaterialTypes = new int[mesh.nElementBoundaries_global];
  mesh.nodeMaterialTypes = new int[mesh.nNodes_global];
  mesh.elementMaterialTypes = new int[mesh.nElements_global];
  for (int i = 0; i < mesh.nElements_global; ++i)
    mesh.elementMaterialTypes[i] = DEFAULT_ELEMENT_MATERIAL;
  int dim = m->getDimension();
  apf::MeshIterator* it = m->begin(dim - 1);
  apf::MeshEntity* f;
  while ((f = m->iterate(it))) {
    int i = localNumber(f);
    mesh.elementBoundaryMaterialTypes[i] = getInOutMaterial(m, f);
  }
  m->end(it);
  it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
    int i = localNumber(v);
    mesh.nodeMaterialTypes[i] = getInOutMaterial(m, v);
  }
  m->end(it);
  return 0;
}

/* Given a geometric model face identified by the integer
 * (scorec_tag), get all nodes and element boundaries
 * classified on the closure of that model face and
 * put the (proteus_material) integer in their material
 * array slot.
 */
int MeshAdaptPUMIDrvr::updateMaterialArrays(Mesh& mesh,
    int dim,
    int proteus_material,
    int scorec_tag)
{
  //int dim = m->getDimension();
  apf::ModelEntity* geomEnt = m->findModelEntity(dim, scorec_tag);
  apf::MeshIterator* it = m->begin(dim);
  apf::MeshEntity* f;
  if(dim==m->getDimension()){
    while ((f = m->iterate(it))) {
     if (m->toModel(f) == geomEnt) {
      int i = localNumber(f);
      mesh.elementMaterialTypes[i] = proteus_material;
     }
    }
  }
  else{
    while ((f = m->iterate(it))) {
      if (m->toModel(f) == geomEnt) {
        int i = localNumber(f);
        mesh.elementBoundaryMaterialTypes[i] = proteus_material;
      }
    }
    apf::DynamicArray<apf::Node> nodes;
    apf::getNodesOnClosure(m, geomEnt, nodes);
    for (size_t i = 0; i < nodes.getSize(); ++i) {
      int vtxId = localNumber(nodes[i].entity);
      mesh.nodeMaterialTypes[vtxId] = proteus_material;
    }
  }
  m->end(it);
  return 0;
}

/* Overload updateMaterialArray for reconstructed SCOREC meshes.
 * The material types are stored based on derived model entity.
 * We can recover the material of a mesh entity by looking at its classified
 * model entity
 */
int MeshAdaptPUMIDrvr::updateMaterialArrays(Mesh& mesh)
{
  std::cout<<"UPDATING MATERIALS FOR RECONSTRUCTED MESH\n";
  int geomTag;
  apf::ModelEntity* geomEnt;
  apf::MeshIterator* it;
  apf::MeshEntity* f;

  int dim = 0;
  it = m->begin(dim);
  while(f = m->iterate(it)){
    int i = localNumber(f);
    geomEnt = m->toModel(f);
    geomTag = m->getModelTag(geomEnt);
    if(m->getModelType(geomEnt) == dim){
      mesh.nodeMaterialTypes[i] =modelVertexMaterial[geomTag];
      std::cout<<"This is the geomTag "<<geomTag<<" this is the material "<<modelVertexMaterial[geomTag]<<std::endl;
    }
    else if(m->getModelType(geomEnt)==(m->getDimension()-1)){ //on the boundary entity
      mesh.nodeMaterialTypes[i] =modelBoundaryMaterial[geomTag];
    }
    else{
      mesh.nodeMaterialTypes[i] = 0; //This assumes that all vertices on the boundary are model vertices
    }
  }
  m->end(it);
  //std::abort();
  if(m->getDimension()==2)
    dim = 1;
  else
    dim = 2;
  it = m->begin(dim);
  while(f = m->iterate(it)){
    geomEnt = m->toModel(f);
    int i = localNumber(f);
    if(m->countUpward(f)==1){//necessarily a boundary mesh edge
      std::cout<<"Edge "<<i<<" geomType "<<m->getModelType(geomEnt)<<" geomTag "<<m->getModelTag(m->toModel(f))<<" material "<<modelBoundaryMaterial[m->getModelTag(m->toModel(f))]<<std::endl;
    }
    if(m->getModelType(geomEnt) == dim){
      geomTag = m->getModelTag(m->toModel(f));
      //std::cout<<"Entity "<<i<<" geomTag "<<geomTag<<" geomType "<<m->getModelType(geomEnt)<<" initial material "<<mesh.elementBoundaryMaterialTypes[i]<<" after "<<modelBoundaryMaterial[geomTag]<<std::endl;
      mesh.elementBoundaryMaterialTypes[i] =modelBoundaryMaterial[geomTag];
    }
    else{
      geomTag = m->getModelTag(m->toModel(f));
      //std::cout<<"Entity "<<i<<" geomTag "<<geomTag<<" geomType "<<m->getModelType(geomEnt)<<" initial material "<<mesh.elementBoundaryMaterialTypes[i]<<" after "<<0<<std::endl;
      

      //THIS LOOKS LIKE A BUG THAT WILL NEED TO BE FIXED. SHOULD BE REGION DEPENDENT
      mesh.elementBoundaryMaterialTypes[i] = 0; 
    }
  }
  m->end(it);

  dim = dim+1; //the regions are necessarily one dimension higher than previous dim
  it = m->begin(dim);
  while(f = m->iterate(it)){
    geomEnt = m->toModel(f);
    int i = localNumber(f);
    assert(m->getModelType(geomEnt)==dim);
    geomTag = m->getModelTag(m->toModel(f));
    //The geomTag is actually the right material for region entities
    mesh.elementMaterialTypes[i] = geomTag; //modelRegionMaterial[geomTag];
  }
  m->end(it);
  return 0;
}


#include <PCU.h>
#include "apfConvert.h"
#include "apfMesh2.h"
#include "apf.h"
#include "apfNumbering.h"
#include <map>

namespace apf {

typedef int Gid;

static void constructVerts(
    Mesh2* m, int nverts,
    int* local2globalMap,
    GlobalToVert& result)
{
  ModelEntity* interior = m->findModelEntity(m->getDimension(), 0);
  for (int i = 0; i < nverts; ++i)
    result[local2globalMap[i]] = m->createVert_(interior);
}


static void constructBoundaryElements(
    Mesh2* m, const Gid* conn_b, int nelem_b, int etype_b,
    GlobalToVert& globalToVert)
{
  ModelEntity* interior = m->findModelEntity(m->getDimension(), 0);
  int nev = apf::Mesh::adjacentCount[etype_b][0];
  for (int i = 0; i < nelem_b; ++i) {
    Downward verts;
    int offset = i * nev;
    for (int j = 0; j < nev; ++j)
      verts[j] = globalToVert[conn_b[j + offset]];
    m->createEntity(etype_b,interior,verts);
  }
}
static void constructElements(
    Mesh2* m, const Gid* conn, int nelem, int etype,
    GlobalToVert& globalToVert)
{
  ModelEntity* interior = m->findModelEntity(m->getDimension(), 0);
  int nev = apf::Mesh::adjacentCount[etype][0];
  for (int i = 0; i < nelem; ++i) {
    Downward verts;
    int offset = i * nev;
    for (int j = 0; j < nev; ++j)
      verts[j] = globalToVert[conn[j + offset]];
    buildElement(m, interior, etype, verts);
  }
}

static Gid getMax(const GlobalToVert& globalToVert)
{
  Gid max = -1;
  APF_CONST_ITERATE(GlobalToVert, globalToVert, it)
    max = std::max(max, it->first);
  return PCU_Max_Int(max); // this is type-dependent
}


/* algorithm courtesy of Sebastian Rettenberger:
   use brokers/routers for the vertex global ids.
   Although we have used this trick before (see mpas/apfMPAS.cc),
   I didn't think to use it here, so credit is given. */
static void constructResidence(Mesh2* m, GlobalToVert& globalToVert)
{
  Gid max = getMax(globalToVert);
  Gid total = max + 1;
  int peers = PCU_Comm_Peers();
  int quotient = total / peers;
  int remainder = total % peers;
  int mySize = quotient;
  int self = PCU_Comm_Self();
  if (self == (peers - 1))
    mySize += remainder;
  typedef std::vector< std::vector<int> > TmpParts;
  TmpParts tmpParts(mySize);
  /* if we have a vertex, send its global id to the
     broker for that global id */
  PCU_Comm_Begin();
  APF_ITERATE(GlobalToVert, globalToVert, it) {
    int gid = it->first;
    int to = std::min(peers - 1, gid / quotient);
    PCU_COMM_PACK(to, gid);
  }
  PCU_Comm_Send();
  int myOffset = self * quotient;
  /* brokers store all the part ids that sent messages
     for each global id */
  while (PCU_Comm_Receive()) {
    int gid;
    PCU_COMM_UNPACK(gid);
    int from = PCU_Comm_Sender();
    tmpParts.at(gid - myOffset).push_back(from);
  }
  /* for each global id, send all associated part ids
     to all associated parts */
  PCU_Comm_Begin();
  for (int i = 0; i < mySize; ++i) {
    std::vector<int>& parts = tmpParts[i];
    for (size_t j = 0; j < parts.size(); ++j) {
      int to = parts[j];
      int gid = i + myOffset;
      int nparts = parts.size();
      PCU_COMM_PACK(to, gid);
      PCU_COMM_PACK(to, nparts);
      for (size_t k = 0; k < parts.size(); ++k)
        PCU_COMM_PACK(to, parts[k]);
    }
  }
  PCU_Comm_Send();
  /* receiving a global id and associated parts,
     lookup the vertex and classify it on the partition
     model entity for that set of parts */
  while (PCU_Comm_Receive()) {
    int gid;
    PCU_COMM_UNPACK(gid);
    int nparts;
    PCU_COMM_UNPACK(nparts);
    Parts residence;
    for (int i = 0; i < nparts; ++i) {
      int part;
      PCU_COMM_UNPACK(part);
      residence.insert(part);
    }
    MeshEntity* vert = globalToVert[gid];
    m->setResidence(vert, residence);
  }
}

/* given correct residence from the above algorithm,
   negotiate remote copies by exchanging (gid,pointer)
   pairs with parts in the residence of the vertex */
static void constructRemotes(Mesh2* m, GlobalToVert& globalToVert)
{
  int self = PCU_Comm_Self();
  PCU_Comm_Begin();
  APF_ITERATE(GlobalToVert, globalToVert, it) {
    int gid = it->first;
    MeshEntity* vert = it->second;
    Parts residence;
    m->getResidence(vert, residence);
    APF_ITERATE(Parts, residence, rit)
      if (*rit != self) {
        PCU_COMM_PACK(*rit, gid);
        PCU_COMM_PACK(*rit, vert);
      }
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    int gid;
    PCU_COMM_UNPACK(gid);
    MeshEntity* remote;
    PCU_COMM_UNPACK(remote);
    int from = PCU_Comm_Sender();
    MeshEntity* vert = globalToVert[gid];
    m->addRemote(vert, from, remote);
  }
}

void construct(Mesh2* m, const int* conn, const int* conn_b, int nelem, 
    int nelem_b, int nverts,int etype, int etype_b, int* local2globalMap,
    GlobalToVert& globalToVert)
{
  std::cout<<"Entering Custom Construct Function\n";
  constructVerts(m, nverts,local2globalMap,globalToVert);
  constructBoundaryElements(m, conn_b, nelem_b, etype_b, globalToVert);
  constructElements(m, conn, nelem, etype, globalToVert);
  constructResidence(m, globalToVert);
  constructRemotes(m, globalToVert);
  PCU_Barrier();
  stitchMesh(m);
  m->acceptChanges();
  std::cout<<"Exiting Custom Construct Function\n";
}

}

//Attempt to reconstruct a PUMI mesh based on Proteus mesh data
////structures.
#include <apf.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <gmi.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfConvert.h>
#include <ma.h>
#include <PCU.h>

#include <cassert>
#include <gmi_lookup.h>

int MeshAdaptPUMIDrvr::transferModelInfo(int* numGeomEntities, int* edges, int* faces, int* mVertex2Model, int*mEdge2Model, int*mBoundary2Model,int nMaxSegments){
  numModelEntities[0] = numGeomEntities[0];
  numModelEntities[1] = numGeomEntities[1];
  numModelEntities[2] = numGeomEntities[2];
  numModelEntities[3] = numGeomEntities[3];
  edgeList = edges;
  faceList = faces;
  meshVertex2Model = mVertex2Model;
  meshEdge2Model = mEdge2Model;
  meshBoundary2Model = mBoundary2Model;
  numSegments = nMaxSegments;

  std::cerr<<"Finished Transferring Model Info\n";
  return 0;
}

int MeshAdaptPUMIDrvr::reconstructFromProteus(Mesh& mesh, Mesh& globalMesh,int hasModel)
{
  std::cout<<"STARTING RECONSTRUCTION!\n";
  isReconstructed = 1; //True

  //Preliminaries
  comm_size = PCU_Comm_Peers();
  comm_rank = PCU_Comm_Self();

  int numModelNodes;
  int numModelEdges;
  int numModelBoundaries;
  int numModelRegions;

  int nBoundaryNodes=0; //number of total boundary nodes regardless of ownership
  int nNodes_owned = globalMesh.nodeOffsets_subdomain_owned[PCU_Comm_Self()+1]-globalMesh.nodeOffsets_subdomain_owned[PCU_Comm_Self()];
  if(PCU_Comm_Self()==0){
    std::cout<<"This is rank "<<0<<" nNodes_owned "<<nNodes_owned<<std::endl;
  }
  PCU_Barrier();
  if(PCU_Comm_Self()==1){
    std::cout<<"This is rank "<<1<<" nNodes_owned "<<nNodes_owned<<std::endl;
  }
  PCU_Barrier();

  for(int i =0;i<mesh.nNodes_global;i++){
    if(mesh.nodeMaterialTypes[i]>0){
      nBoundaryNodes++;
    }    
  }

  int numDim;
  if(mesh.nNodes_element==3)
    numDim = 2;
  else
    numDim = 3;
  if(hasModel){
    numModelNodes=numModelEntities[0];
    numModelEdges=numModelEntities[1];
    numModelBoundaries=numModelEntities[2];
    numModelRegions=numModelEntities[3];
    if(numDim=2){
      //should add some sort of assertion statement here
      numModelBoundaries = numModelEdges;
    }
  }
  else{
    numModelNodes = nBoundaryNodes;
    numModelEdges = mesh.nEdges_global;
    numModelBoundaries = mesh.nExteriorElementBoundaries_global;
    numModelRegions = numModelEntities[3];
  }

  assert(numModelRegions>0);
  //create Model
  gmi_model* gMod;

  struct gmi_base* gMod_base;
  gMod_base = (gmi_base*)malloc(sizeof(*gMod_base));
  gMod_base->model.ops = &gmi_base_ops;
  gmi_base_init(gMod_base);

  struct agm_ent e;
  struct agm_bdry b;
  struct agm_ent d;

  PCU_Barrier();
  if(comm_rank==0){
    std::cout<<"This is comm rank BEFORE SYNC "<<PCU_Comm_Self()<<std::endl;
    std::cout<<"Number of model entities: "<<numModelNodes<<" "<<numModelEdges<<" "<<numModelBoundaries<<" "<<numModelRegions<<std::endl;
    //std::cout<<"Number of mesh entities: "<<mesh.nNodes_owned << " "<<mesh.nEdges_owned<<" "<<mesh.nElementBoundaries_owned<<" "<<mesh.nElements_owned<<std::endl;
    std::cout<<"Number of mesh entities: "<<mesh.nNodes_global << " "<<mesh.nEdges_global<<" "<<mesh.nElementBoundaries_global<<" "<<mesh.nElements_global<<std::endl;
  }
  PCU_Barrier();
  if(comm_rank==1){
    std::cout<<"This is comm rank BEFORE SYNC "<<PCU_Comm_Self()<<std::endl;
    std::cout<<"Number of model entities: "<<numModelNodes<<" "<<numModelEdges<<" "<<numModelBoundaries<<" "<<numModelRegions<<std::endl;
    //std::cout<<"Number of mesh entities: "<<mesh.nNodes_owned << " "<<mesh.nEdges_owned<<" "<<mesh.nElementBoundaries_owned<<" "<<mesh.nElements_owned<<std::endl;
    std::cout<<"Number of mesh entities: "<<mesh.nNodes_global << " "<<mesh.nEdges_global<<" "<<mesh.nElementBoundaries_global<<" "<<mesh.nElements_global<<std::endl;
  }
  PCU_Barrier();

  numModelTotals[0] = numModelNodes;
  numModelTotals[1] = numModelEdges;
  numModelTotals[2] = numModelBoundaries;
  numModelTotals[3] = 0;//The total number of regions is known so no need to set a value
  PCU_Add_Ints(&numModelTotals[0],4); //get all offsets at the same time
  numModelTotals[3] = numModelRegions;

  //gvertices
  gmi_base_reserve(gMod_base,AGM_VERTEX,numModelTotals[0]);

  if(numDim==2){
  //gedges
    gmi_base_reserve(gMod_base,AGM_EDGE,numModelTotals[2]);
  //gfaces
    gmi_base_reserve(gMod_base,AGM_FACE,numModelTotals[3]);
  //gregions
    gmi_base_reserve(gMod_base,AGM_REGION,0);
  }
  else if(numDim==3){
  //gedges
    gmi_base_reserve(gMod_base,AGM_EDGE,numModelTotals[1]);
  //gfaces
    gmi_base_reserve(gMod_base,AGM_FACE,numModelTotals[2]);
  //gregions
    gmi_base_reserve(gMod_base,AGM_REGION,numModelTotals[3]);
  }

  gMod = &gMod_base->model;

  std::cout<<"Passed the model allocation\n";
  //create Mesh
  m = apf::makeEmptyMdsMesh(gMod,2,false);

  //We can use apf::construct() and a set coordinates function to generate the mesh
  
  int etype,etype_b;
  apf::GlobalToVert outMap;
  if(numDim == 2){
    etype = apf::Mesh::TRIANGLE;
    etype_b = apf::Mesh::EDGE;
  }
  else{
    etype = apf::Mesh::TET;
    etype_b = apf::Mesh::TRIANGLE;
  }

  std::cout<<"At construction site\n";

  int* local2global_elementBoundaryNodes;
  local2global_elementBoundaryNodes = (int*) malloc(sizeof(int)*mesh.nElementBoundaries_global*apf::Mesh::adjacentCount[etype_b][0]);
  for(int i=0;i<mesh.nElementBoundaries_global*apf::Mesh::adjacentCount[etype_b][0];i++){ //should use adjacent count function from core
    local2global_elementBoundaryNodes[i] = globalMesh.nodeNumbering_subdomain2global[mesh.elementBoundaryNodesArray[i]];
  }
  PCU_Barrier();
  int* local2global_elementNodes;
  local2global_elementNodes = (int*) malloc(sizeof(int)*mesh.nElements_global*apf::Mesh::adjacentCount[etype][0]);
  for(int i=0;i<mesh.nElements_global*apf::Mesh::adjacentCount[etype][0];i++){ //should use adjacent count function from core
    local2global_elementNodes[i] = globalMesh.nodeNumbering_subdomain2global[mesh.elementNodesArray[i]];
  }


  apf::construct(m,local2global_elementNodes,local2global_elementBoundaryNodes,
    mesh.nElements_global,mesh.nElementBoundaries_global,mesh.nNodes_global,etype,etype_b,
    globalMesh.nodeNumbering_subdomain2global,outMap);

  PCU_Barrier();
  std::cout<<"Past construction site\n";

  //Get the global model offsets after the mesh has been created
  //Need to get the number of owned element boundaries on the current rank
  //Also need to get the number of owned exterior entities for proper processor communication

  nBoundaryNodes = 0;
  apf::MeshIterator* entIter=m->begin(0);
  apf::MeshEntity* ent;
  int idx = 0;
  while((ent=m->iterate(entIter))){
    if(m->isOwned(ent)){
      if(mesh.nodeMaterialTypes[idx]>0)
        nBoundaryNodes++;
    }
    idx++;
  }
  m->end(entIter);

  //int nElementBoundaries_owned = globalMesh.elementBoundaryOffsets_subdomain_owned[PCU_Comm_Self()+1]-globalMesh.elementBoundaryOffsets_subdomain_owned[PCU_Comm_Self()];
  entIter=m->begin(1);
  idx=0;
  int nExteriorElementBoundaries_owned = 0;
  while((ent=m->iterate(entIter))){
    if(m->isOwned(ent)){
      if(mesh.elementBoundaryMaterialTypes[idx]>0)
        nExteriorElementBoundaries_owned++;  
    }
    idx++;
  }
  m->end(entIter);

  if(hasModel){
    numModelNodes=numModelEntities[0];
    numModelEdges=numModelEntities[1];
    numModelBoundaries=numModelEntities[2];
    numModelRegions=numModelEntities[3];
    if(numDim=2){
      //should add some sort of assertion statement here
      numModelBoundaries = numModelEdges;
    }
  }
  else{
    numModelNodes = nBoundaryNodes;
    numModelEdges = mesh.nEdges_global;
    numModelBoundaries = nExteriorElementBoundaries_owned;
    numModelRegions = numModelEntities[3]; 
  }

  //////////

  numModelOffsets[0] = numModelNodes;
  numModelOffsets[1] = numModelEdges;
  numModelOffsets[2] = numModelBoundaries;
  numModelOffsets[3] = 0;//numModelRegions; what happens with multiple regions?
  
  numModelTotals[0] = numModelNodes;
  numModelTotals[1] = numModelEdges;
  numModelTotals[2] = numModelBoundaries;
  numModelTotals[3] = 0;//numModelRegions; what happens with multiple regions?

  PCU_Barrier();
  PCU_Exscan_Ints(&numModelOffsets[0],4); //get all offsets at the same time
  PCU_Add_Ints(&numModelTotals[0],4); //get all offsets at the same time
  numModelTotals[3] = numModelRegions;
  PCU_Barrier();
  if(comm_rank==0){
    for(int i=0;i<4;i++){
      std::cout<<"What is the model offsets for this rank? "<<numModelOffsets[i]<<std::endl;
      std::cout<<"What is the model total for this rank? "<<numModelTotals[i]<<std::endl;
    }
  }
  PCU_Barrier();
  if(comm_rank==1){
    for(int i=0;i<4;i++){
      std::cout<<"What is the model offsets for this rank? "<<numModelOffsets[i]<<std::endl;
      std::cout<<"What is the model total for this rank? "<<numModelTotals[i]<<std::endl;
    }
  }
  PCU_Barrier();
  if(comm_rank==0){
    std::cout<<"This is comm rank "<<PCU_Comm_Self()<<std::endl;
    std::cout<<"Number of model entities: "<<numModelNodes<<" "<<numModelEdges<<" "<<numModelBoundaries<<" "<<numModelRegions<<std::endl;
    //std::cout<<"Number of mesh entities: "<<mesh.nNodes_owned << " "<<mesh.nEdges_owned<<" "<<mesh.nElementBoundaries_owned<<" "<<mesh.nElements_owned<<std::endl;
    std::cout<<"Number of mesh entities: "<<mesh.nNodes_global << " "<<mesh.nEdges_global<<" "<<mesh.nElementBoundaries_global<<" "<<mesh.nElements_global<<std::endl;
  }
  PCU_Barrier();
  if(comm_rank==1){
    std::cout<<"This is comm rank "<<PCU_Comm_Self()<<std::endl;
    std::cout<<"Number of model entities: "<<numModelNodes<<" "<<numModelEdges<<" "<<numModelBoundaries<<" "<<numModelRegions<<std::endl;
    //std::cout<<"Number of mesh entities: "<<mesh.nNodes_owned << " "<<mesh.nEdges_owned<<" "<<mesh.nElementBoundaries_owned<<" "<<mesh.nElements_owned<<std::endl;
    std::cout<<"Number of mesh entities: "<<mesh.nNodes_global << " "<<mesh.nEdges_global<<" "<<mesh.nElementBoundaries_global<<" "<<mesh.nElements_global<<std::endl;
  }
  PCU_Barrier();

  //classify mesh entities on model entities

  apf::Vector3 pt;

  apf::ModelEntity* g_vertEnt;
  apf::ModelEntity* g_edgeEnt;
  apf::ModelEntity* g_faceEnt;
  apf::MeshEntity* vertEnt;
  
  modelVertexMaterial = (int*)calloc(numModelTotals[0],sizeof(int));
  modelBoundaryMaterial = (int*)calloc(numModelTotals[2],sizeof(int));
  modelRegionMaterial = (int*)calloc(numModelTotals[3],sizeof(int));

  //gmi set entities

  gmi_unfreeze_lookups(gMod_base->lookup);
  for(int i=0;i<numModelTotals[0];i++){
    e = agm_add_ent(gMod_base->topo, AGM_VERTEX);
    gmi_set_lookup(gMod_base->lookup, e, i);
  }
  gmi_freeze_lookup(gMod_base->lookup, (agm_ent_type)0);

  for(int i=0;i<numModelTotals[2];i++){
    e = agm_add_ent(gMod_base->topo, AGM_EDGE);
    gmi_set_lookup(gMod_base->lookup, e, i);
  }
  gmi_freeze_lookup(gMod_base->lookup, (agm_ent_type)1);

  for(int i=0;i<numModelRegions;i++){
    e = agm_add_ent(gMod_base->topo, AGM_FACE);
    gmi_set_lookup(gMod_base->lookup, e, i+1); //assumes material types are enumerated starting from 1
    if(hasModel){
      b = agm_add_bdry(gMod_base->topo, e);
      for(int k=0;k<numSegments;k++){
        if(faceList[(i)*numSegments+k]==-1) break;
        else{
          std::cout<<"edge "<<faceList[(i)*numSegments+k]<<" "<<numSegments<<std::endl;
          d = gmi_look_up(gMod_base->lookup,AGM_EDGE,faceList[(i)*numSegments+k]);
          agm_add_use(gMod_base->topo,b,d);
        }
      }
    }
  }
  gmi_freeze_lookup(gMod_base->lookup, (agm_ent_type)2);
  std::cout<<"Finished creating model entities\n";
  PCU_Barrier();


  int matTag; //mesh.elementMaterialType[fID];
  apf::ModelEntity* gEnt; 
  int vertCounter = numModelOffsets[0];//0;

  //Iterate over the vertices and set the coordinates if owned
  int vID = 0;
  entIter = m->begin(0);
  PCU_Comm_Begin();
  while(ent = m->iterate(entIter)){
    //int vID = it->first;
    pt[0]=mesh.nodeArray[vID*3+0];
    pt[1]=mesh.nodeArray[vID*3+1];
    pt[2]=mesh.nodeArray[vID*3+2];
    m->setPoint(ent,0,pt);
    if(m->isOwned(ent)){
      matTag = mesh.nodeMaterialTypes[vID];
      std::cout<<"What is the material type? "<<matTag<<std::endl;
      if(hasModel){
        gEnt = m->findModelEntity(meshVertex2Model[2*vID+1],meshVertex2Model[2*vID]);
        if(meshVertex2Model[2*vID+1]==0) //if entity is a model vertex
          modelVertexMaterial[meshVertex2Model[2*vID]] = matTag;
      }
      else{
        if(matTag==0){
          matTag = mesh.elementMaterialTypes[mesh.nodeElementsArray[mesh.nodeElementOffsets[vID]]];
          gEnt = m->findModelEntity(2,matTag);
        }
        else{
          gEnt = m->findModelEntity(0,vertCounter);
          modelVertexMaterial[vertCounter] = matTag;
          vertCounter++;  
        }
      }
      m->setModelEntity(ent,gEnt);
      if(m->isShared(ent)){
        apf::Copies remotes;
        m->getRemotes(ent,remotes);
        for(apf::Copies::iterator it = remotes.begin(); it != remotes.end(); ++it){
          PCU_COMM_PACK(it->first,it->second);
          PCU_COMM_PACK(it->first,gEnt);
        }
      }
    } //endif owned
    vID++;
  }
  PCU_Comm_Send();
  //receive model entity classification from owning nodes
  while(PCU_Comm_Receive()){
    PCU_COMM_UNPACK(ent);
    PCU_COMM_UNPACK(gEnt);
    m->setModelEntity(ent,gEnt); 
  }
  PCU_Barrier();
  m->end(entIter);

  std::cout<<"Finished setting entities "<<PCU_Comm_Self()<<std::endl;
  PCU_Barrier();

  apf::writeVtkFiles("initialConstructedMesh", m);

  PCU_Barrier();
  if(PCU_Comm_Self()==0){
    std::cout<<"How many verts? "<<m->count(0)<<" owned? "<<apf::countOwned(m,0)<<std::endl;
    std::cout<<"what is vertcounter? "<<vertCounter<<std::endl;
    std::cout<<"what is numModelNodes? "<<numModelNodes<<std::endl;
  }
  PCU_Barrier();
  if(PCU_Comm_Self()==1){
    std::cout<<"what is vertcounter? "<<vertCounter<<std::endl;
    std::cout<<"what is numModelNodes? "<<numModelNodes<<std::endl;
  }
  PCU_Barrier();

  //Classify the mesh edge entities
  //If the edge is on a model edge, it should have a material tag greater than 0.
  //If the edge is on a partition boundary, the material tag should be 0.
  //There is no direct control over ownership when constructing the mesh, so it
  //must be left general.
  int edgID = 0;
  int edgCounter = 0;
  int edgMaterialCounter = numModelOffsets[2];
  entIter=m->begin(1);
  while(ent = m->iterate(entIter)){
    if(PCU_Comm_Self()==1){
      std::cout<<"edge counter "<<edgCounter<<" boundary "<< mesh.exteriorElementBoundariesArray[edgCounter]<<" edgID "<<edgID<<" material "<<mesh.elementBoundaryMaterialTypes[edgID]<<std::endl;
    }
    if(hasModel){
      gEnt = m->findModelEntity(meshBoundary2Model[2*edgID+1],meshBoundary2Model[2*edgID]);
      //std::cout<<"What is the search say? "<<m->getModelType(gEnt)<<" "<<m->getModelTag(gEnt)<<" Type "<<meshBoundary2Model[2*edgID+1]<<" ID "<<meshBoundary2Model[2*edgID]<<std::endl;
      if(meshBoundary2Model[2*edgID+1]==1) //if entity is a on a model boundary
        modelBoundaryMaterial[meshBoundary2Model[2*edgID]] = mesh.elementBoundaryMaterialTypes[edgID];
    }
    else{
      if(mesh.exteriorElementBoundariesArray[edgCounter]==edgID && (mesh.elementBoundaryMaterialTypes[edgID]!=0)){
        gEnt = m->findModelEntity(1,edgMaterialCounter);
        modelBoundaryMaterial[edgMaterialCounter] = mesh.elementBoundaryMaterialTypes[edgID]; 
        edgCounter++;
        edgMaterialCounter++;
      }
      else {
        //If the current exterior entity is an edge on a partition boundary, need to check material and
        //get to the next item in the exterior array
        if(m->isShared(ent) && mesh.elementBoundaryMaterialTypes[mesh.exteriorElementBoundariesArray[edgCounter]]==0) 
          edgCounter++; 
        if(mesh.elementBoundaryMaterialTypes[edgID]!=0){
          if(PCU_Comm_Self()==1){
          std::cout<<edgCounter<<" Failing rank is "<<PCU_Comm_Self()<<" "<<edgID<<" coordinates "<<mesh.elementBoundaryBarycentersArray[edgID*3+0]<<" "<<
            mesh.elementBoundaryBarycentersArray[edgID*3+1]<<" material "<<mesh.elementBoundaryMaterialTypes[edgID]<<
            " "<<local2global_elementBoundaryNodes[2*edgID]<<" "<<local2global_elementBoundaryNodes[2*edgID+1]<<std::endl;
          }
        }
        //There are always two entities adjacent to an element boundary
        //Pick one and take that as the material type for classification
        matTag = mesh.elementMaterialTypes[mesh.elementBoundaryElementsArray[2*edgID]];
        gEnt = m->findModelEntity(2,matTag);
      }
    }
    m->setModelEntity(ent,gEnt);
    edgID++;
  }
  m->end(entIter);

  entIter = m->begin(2);
  PCU_Barrier();
  std::cout<<"Initializing RECONSTRUCTION!\n";

  //Populate the region materials
  //Assumes that the regions are numbered sequentially from 1 onward
  for(int i=0;i<numModelRegions;i++)
    modelRegionMaterial[i] = i+1;
  
  int fID = 0;
  while(ent = m->iterate(entIter)){
    gEnt = m->findModelEntity(2,mesh.elementMaterialTypes[fID]);
    m->setModelEntity(ent,gEnt);
    if(fID==215){
      std::cout<<"THIS IS THE MATERIAL for 215 "<<mesh.elementMaterialTypes[fID]<<std::endl;
      std::cout<<"This is the model entity "<<m->getModelTag(m->toModel(ent))<<std::endl;
    }

    fID++;
  }
  m->end(entIter);

  //Still need to add all of the arrays together for total model material arrays
  
  PCU_Add_Ints(modelVertexMaterial,numModelTotals[0]);
  PCU_Add_Ints(modelBoundaryMaterial,numModelTotals[2]);
  PCU_Add_Ints(modelRegionMaterial,numModelTotals[3]);

  free(local2global_elementBoundaryNodes);
  free(local2global_elementNodes);
  m->acceptChanges();
  apf::alignMdsRemotes(m);
  m->verify();
  initialReconstructed = 1;
  numberLocally();
  m->verify();
  apf::writeVtkFiles("reconstructedMesh", m);
  m->writeNative("constructedMesh_withAdaptIssues.smb");
  gmi_write_dmg(gMod,"constructedMesh_withAdaptIssues.dmg");
  std::cout<<"COMPLETED RECONSTRUCTION!\n";
}



