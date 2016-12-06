#include <algorithm>

#include "MeshAdaptPUMI.h"
#include "mesh.h"
#include <apfShape.h>

#include <sstream>

static apf::Numbering* numberOwnedEntitiesFirst(apf::Mesh* m, int dimension)
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
  while ((e = m->iterate(it)))
    if (m->isOwned(e))
      apf::number(n, e, 0, 0, i++);
  m->end(it);
  it = m->begin(dimension);
  while ((e = m->iterate(it)))
    if (!m->isOwned(e))
      apf::number(n, e, 0, 0, i++);
  m->end(it);
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
    local[d] = numberOwnedEntitiesFirst(m, d);
  }
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
    int proteus_material,
    int scorec_tag)
{
  int dim = m->getDimension();
  apf::ModelEntity* geomEnt = m->findModelEntity(dim - 1, scorec_tag);
  apf::MeshIterator* it = m->begin(dim - 1);
  apf::MeshEntity* f;
  while ((f = m->iterate(it))) {
    if (m->toModel(f) == geomEnt) {
      int i = localNumber(f);
      mesh.elementBoundaryMaterialTypes[i] = proteus_material;
    }
  }
  m->end(it);
  apf::DynamicArray<apf::Node> nodes;
  apf::getNodesOnClosure(m, geomEnt, nodes);
  for (size_t i = 0; i < nodes.getSize(); ++i) {
    int vtxId = localNumber(nodes[i].entity);
    mesh.nodeMaterialTypes[vtxId] = proteus_material;
  }
  return 0;
}
