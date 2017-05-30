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
  int geomTag;
  apf::ModelEntity* geomEnt;
  apf::MeshIterator* it;
  apf::MeshEntity* f;

  int dim = 0;
  it = m->begin(dim);
  while(f = m->iterate(it)){
    geomEnt = m->toModel(f);
    int i = localNumber(f);
    geomTag = m->getModelTag(m->toModel(f));
    if(m->getModelType(geomEnt) == dim){
      mesh.nodeMaterialTypes[i] =modelVertexMaterial[geomTag];
    }
    else if(m->getModelType(geomEnt)==(m->getDimension()-1)){ //on the boundary entity
      mesh.nodeMaterialTypes[i] =modelBoundaryMaterial[geomTag];
    }
    else{
      mesh.nodeMaterialTypes[i] = 0; //This assumes that all vertices on the boundary are model vertices
    }
  }
  m->end(it);

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
      mesh.elementBoundaryMaterialTypes[i] = 0; //This assumes that all vertices on the boundary are model vertices
    }
  }
  m->end(it);

  dim = dim+1; //the regions are necessary one dimension higher than previous dim
  it = m->begin(dim);
  while(f = m->iterate(it)){
    geomEnt = m->toModel(f);
    int i = localNumber(f);
    assert(m->getModelType(geomEnt)==dim);
    geomTag = m->getModelTag(m->toModel(f));
    mesh.elementMaterialTypes[i] =modelRegionMaterial[0];
  }
  m->end(it);

  return 0;
}



//Attempt to reconstruct a PUMI mesh based on Proteus mesh data
////structures.
#include <apf.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <gmi.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <ma.h>

#include <cassert>
#include <gmi_lookup.h>

#include<iostream>
//#include<CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Gmpq.h>
#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Point_3.h>
//#include <CGAL/Segment_3.h>

//typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Simple_cartesian<CGAL::Gmpq> Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Segment_3 Segment_3;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Segment_2 Segment_2;

int MeshAdaptPUMIDrvr::getMesh2ModelClassification(Mesh& mesh){

  //initialize checklists
  std::cout<<"mesh.nNodes_global "<<mesh.nNodes_global<<std::endl;
  std::cout<<"mesh.nNodes_element "<<mesh.nNodes_element<<std::endl;
  std::cout<<"mesh.nNodes_elementBoundary "<<mesh.nNodes_elementBoundary<<std::endl;
  std::cout<<"mesh.nElements_global "<<mesh.nElements_global<<std::endl;
  std::cout<<"mesh.nElementBoundaries_global "<<mesh.nElementBoundaries_global<<std::endl;
  checklist_b = (int*) malloc(sizeof(int)*mesh.nElementBoundaries_global);
  meshBoundary2Model = (int*) malloc(sizeof(int)*mesh.nElementBoundaries_global*2);
  
  //Construct the model vertices
  //std::vector<Point_3> modelPoints;
  std::vector<Point_2> modelPoints;
  for(int i=0; i<numModelEntities[0];i++){
    //modelPoints.push_back(Point_3(vertexList[3*i+0],vertexList[3*i+1],vertexList[3*i+2]));
    modelPoints.push_back(Point_2(vertexList[3*i+0],vertexList[3*i+1]));
  } 
  //Construct the model boundaries
  //std::vector<Segment_3> modelBoundaries;
  std::vector<Segment_2> modelBoundaries;
  int numDim = 2;
  int numBoundaries = numModelEntities[numDim-1];
  for(int i=0; i<numBoundaries;i++){
    std::cout<<i<<" model points "<<edgeList[numDim*i+0]<<" "<<edgeList[numDim*i+1]<<std::endl;
    //modelBoundaries.push_back(Segment_3(modelPoints[edgeList[numDim*i+0]],modelPoints[edgeList[numDim*i+1]]));
    modelBoundaries.push_back(Segment_2(modelPoints[edgeList[numDim*i+0]],modelPoints[edgeList[numDim*i+1]]));
  }
  
  //Construct a more exact barycenters array
  //std::vector<Point_3> elementBoundaryCentroidsArray;
  //std::vector<Point_3> centroid_temp;
  std::vector<Point_2> elementBoundaryCentroidsArray;
  std::vector<Point_2> centroid_temp;

  for(int ebN=0;ebN<mesh.nElementBoundaries_global;ebN++){
    for(int nN=0; nN<mesh.nNodes_elementBoundary;nN++){
      int idx = mesh.elementBoundaryNodesArray[mesh.nNodes_elementBoundary*ebN+nN];
      //std::cout<<ebN<<" idx "<<idx<<" "<<mesh.nodeArray[3*idx+0]<<" "<<mesh.nodeArray[3*idx+1]<<mesh.nodeArray[3*idx+2]<<std::endl;
      //centroid_temp.push_back(Point_3(mesh.nodeArray[3*idx+0],mesh.nodeArray[3*idx+1],mesh.nodeArray[3*idx+2]));
      centroid_temp.push_back(Point_2(mesh.nodeArray[3*idx+0],mesh.nodeArray[3*idx+1]));
    }
    //Point_3 centroid;
    Point_2 centroid;
    if(numDim==2){
      centroid = CGAL::midpoint(centroid_temp[0],centroid_temp[1]);
    }
    else{ //3D
      std::abort();
    }
    assert(centroid_temp.size()==mesh.nNodes_elementBoundary);
    centroid_temp.clear();
    //std::cout<<"What is the mesh boundary barycenter "<<barycenter<<std::endl;
    elementBoundaryCentroidsArray.push_back(centroid);
  }
  std::cout<<"End construction of barycenters\n";

  //Loop through external element boundaries
  for(int i=0;i<mesh.nExteriorElementBoundaries_global;i++){
    int idx = mesh.exteriorElementBoundariesArray[i];
    //Point_3 testPoint = Point_3(mesh.elementBoundaryBarycentersArray[3*idx+0],mesh.elementBoundaryBarycentersArray[3*idx+1],mesh.elementBoundaryBarycentersArray[3*idx+2]);
    //Point_3 testPoint = elementBoundaryCentroidsArray[idx];
    Point_2 testPoint = elementBoundaryCentroidsArray[idx];
    for(int idxBoundary=0;idxBoundary<numBoundaries;idxBoundary++){
      //Segment_3 modelBoundary = modelBoundaries[idxBoundary]; 
      Segment_2 modelBoundary = modelBoundaries[idxBoundary]; 
      if(idx == 0){
        std::cout<<idx<<" What is modelBoundary ("<<modelBoundary.vertex(0).x()<<","<<modelBoundary.vertex(1).y()<<"),("<<modelBoundary.vertex(1)<<"),("<<testPoint<<")"<<std::endl;
        std::cout<<"What is the orientation? "<<CGAL::orientation(modelBoundary.vertex(0),testPoint,modelBoundary.vertex(1))<<std::endl;
/*
        std::cout<<"is on? "<<modelBoundary.has_on(testPoint)<<std::endl;
        Point_3 testPoint2 = Point_3(3.137,1.566,0);
        std::cout<<"is on2? "<<modelBoundary.has_on(testPoint2)<<std::endl;
        Point_3 testPoint3 = CGAL::midpoint(modelBoundary.vertex(0),modelBoundary.vertex(1));
        std::cout<<"is on3? "<<modelBoundary.has_on(testPoint3)<<" "<<testPoint3<<std::endl;  
*/
      }
      if(modelBoundary.has_on(testPoint)){
        meshBoundary2Model[2*idx+0] = idxBoundary;
        meshBoundary2Model[2*idx+1] = numDim-1;
        checklist_b[idx] = 1; 
        for(int j=0; j<numDim;j++){
          int vID=mesh.elementBoundaryNodesArray[numDim*idx+j];
          if(!checklist_v[vID]){
            meshVertex2Model[2*vID+0] = idxBoundary;
            meshVertex2Model[2*vID+1] = numDim-1;
            checklist_v[vID] = 1;
          }
        }
      } //end if has_on
    }
  }
  for(int i=0;i<mesh.nInteriorElementBoundaries_global;i++){
    int idx = mesh.interiorElementBoundariesArray[i];
    meshBoundary2Model[2*idx+0] = mesh.elementMaterialTypes[0];   
    meshBoundary2Model[2*idx+1] = numDim;
    checklist_b[idx] = 1;
    for(int j=0; j<numDim;j++){
      int vID=mesh.elementBoundaryNodesArray[numDim*idx+j];
      if(!checklist_v[vID]){
        meshVertex2Model[2*vID+0] = mesh.elementMaterialTypes[0];
        meshVertex2Model[2*vID+1] = numDim;
        checklist_v[vID] = 1;
      }
    }
  }
  for(int i=0;i<mesh.nInteriorElementBoundaries_global;i++){
    if(checklist_b[i]!=1){
      std::cout<<"checklist_b "<<i<<" "<<checklist_b[i]<<std::endl;
    }
  }
  //free(modelPoints);
  //free(modelBoundaries);
  //free(meshBoundary2Model);
  //free(checklist_b);
  std::cerr<<"Reached the end of mesh2model classification\n";
}

int MeshAdaptPUMIDrvr::transferModelInfo(int* numGeomEntities, double* vertices, int* edges, int* faces,int* mVertex2Model, int* vertexChecklist){
  numModelEntities[0] = numGeomEntities[0];
  numModelEntities[1] = numGeomEntities[1];
  numModelEntities[2] = numGeomEntities[2];
  numModelEntities[3] = numGeomEntities[3];
  std::cout<<"How many v e b r "<< numModelEntities[0]<<" "<<numModelEntities[1]<<" "<<numModelEntities[2]<<" "<<numModelEntities[3]<<std::endl;

  vertexList = vertices;
  edgeList = edges;
  faceList = faces;
  meshVertex2Model = mVertex2Model;
  checklist_v = vertexChecklist;
  for(int i=0;i<numModelEntities[0];i++){
    std::cout<<"vertex "<<i<<" "<<vertexList[3*i+0]<<" "<<vertexList[3*i+1]<<std::endl;
  }
  std::cout<<"What is meshVertex2Model? "<<meshVertex2Model[0]<<" "<<meshVertex2Model[1]<<" "<<meshVertex2Model[2]<<std::endl;
  std::cerr<<"Finished Transferring Model Info\n";
  return 0;
}

int MeshAdaptPUMIDrvr::reconstructFromProteus(Mesh& mesh, int hasModel)
{
  std::cout<<"STARTING RECONSTRUCTION!\n";
  isReconstructed = 1; //True

  //Preliminaries
  int numModelNodes;
  int numModelEdges;
  int numModelBoundaries;
  int numModelRegions;

  int nBoundaryNodes=0;
  for(int i =0;i<mesh.nNodes_global;i++){
    if(mesh.nodeMaterialTypes[i]>0){
      nBoundaryNodes++;
    }    
  }

  if(hasModel){
    numModelNodes=numModelEntities[0];
    numModelEdges=numModelEntities[1];
    numModelBoundaries=numModelEntities[2];
    numModelRegions=numModelEntities[3];
    if(numModelBoundaries==0){
      //should add some sort of assertion statement here
      numModelBoundaries = numModelEdges;
    }
  }
  else{
    numModelNodes = nBoundaryNodes;
    numModelEdges = mesh.nEdges_global;
    numModelBoundaries = mesh.nExteriorElementBoundaries_global;
    numModelRegions = 1;
  }
    
  std::cout<<"Number of model entities: "<<numModelNodes<<" "<<numModelEdges<<" "<<numModelBoundaries<<" "<<numModelRegions<<std::endl;
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
  
  //gvertices
  gmi_base_reserve(gMod_base,AGM_VERTEX,numModelNodes);
  for(int i=0;i<numModelNodes;i++){
    e = agm_add_ent(gMod_base->topo, AGM_VERTEX);
    gmi_set_lookup(gMod_base->lookup, e, i);
  }
  gmi_freeze_lookup(gMod_base->lookup, (agm_ent_type)0);

  //gedges
  gmi_base_reserve(gMod_base,AGM_EDGE,numModelEdges);
  for(int i=0;i<numModelEdges;i++){
    e = agm_add_ent(gMod_base->topo, AGM_EDGE);
    gmi_set_lookup(gMod_base->lookup, e, i);
  }
  gmi_freeze_lookup(gMod_base->lookup, (agm_ent_type)1);

  //gfaces
  gmi_base_reserve(gMod_base,AGM_FACE,numModelRegions);

  e = agm_add_ent(gMod_base->topo, AGM_FACE);
  gmi_set_lookup(gMod_base->lookup, e, mesh.elementMaterialTypes[0]); //assumes material types are uniform

  gmi_freeze_lookup(gMod_base->lookup, (agm_ent_type)2);

  gMod = &gMod_base->model;

  //create Mesh
  m = apf::makeEmptyMdsMesh(gMod,2,false);
  apf::Vector3 pt;
  apf::MeshEntity* vertices[mesh.nNodes_global]; //possibly make this a growing list..
  apf::MeshEntity* elemverts[mesh.nNodes_element];
  apf::MeshEntity* edges[mesh.nEdges_global];
  apf::MeshEntity* elements[mesh.nElements_global];

  apf::ModelEntity* g_vertEnt;
  apf::ModelEntity* g_edgeEnt;
  apf::ModelEntity* g_faceEnt;
  apf::MeshEntity* vertEnt;
  
  //these data structures assume that every boundary entity is a model entity
  modelVertexMaterial = (int*)malloc(numModelNodes*sizeof(int));
  modelBoundaryMaterial = (int*)malloc(numModelBoundaries*sizeof(int));
  modelRegionMaterial = (int*)malloc(numModelRegions*sizeof(int));
    
  int matTag; //mesh.elementMaterialType[fID];
  apf::ModelEntity* gEnt; 
  for(int vID=0; vID<mesh.nNodes_global;vID++){
    matTag = mesh.nodeMaterialTypes[vID];
    if(hasModel){
      gEnt = m->findModelEntity(meshVertex2Model[2*vID+1],meshVertex2Model[2*vID]);
      if(meshVertex2Model[2*vID+1]==0) //if entity is a model vertex
        modelVertexMaterial[meshVertex2Model[2*vID]] = matTag;
    }
    else{
      int vertCounter = 0;
      if(matTag==0){
        matTag = mesh.elementMaterialTypes[0]; //also assumes that there is only a single material type
        gEnt = m->findModelEntity(2,matTag);
      }
      else{
        gEnt = m->findModelEntity(0,vertCounter);
        modelVertexMaterial[vertCounter] = matTag;
        vertCounter++;  
      }
    }

    pt[0]=mesh.nodeArray[vID*3+0];
    pt[1]=mesh.nodeArray[vID*3+1];
    pt[2]=mesh.nodeArray[vID*3+2];
    //std::cout<<"Mesh Vertex Dimension is "<<m->getModelType(gEnt)<<std::endl;
    vertices[vID] = m->createVert(gEnt);
    m->setPoint(vertices[vID],0,pt);
  }

  for(int i =0;i<numModelNodes;i++){
    std::cout<<"Model vertex "<< i << " has tag "<< modelVertexMaterial[i]<<std::endl;
  }
    
  //Construct the mesh edge entities
  
  apf::Downward down_edge;
/*
  for(int edgID=0; edgID<mesh.nExteriorElementBoundaries_global;edgID++){
    gEnt = m->findModelEntity(1,edgID);
    modelBoundaryMaterial[edgID] = mesh.elementBoundaryMaterialTypes[edgID]; 
    std::cout<<edgID<<" What is the boundary material? "<<mesh.elementBoundaryMaterialTypes[edgID]<<std::endl;
    //i need the edge adjacency array to tell me what the adjacent edges should be..
    int actual_edgID = mesh.exteriorElementBoundariesArray[edgID];
    down_edge[0] = vertices[mesh.edgeNodesArray[actual_edgID*2]];
    down_edge[1] = vertices[mesh.edgeNodesArray[actual_edgID*2+1]];
    //std::cout<<edgID<<" Adjacent vertices are "<< mesh.edgeNodesArray[actual_edgID*2]<<" "<<mesh.edgeNodesArray[actual_edgID*2+1]<<std::endl;
    m->createEntity(apf::Mesh::EDGE,gEnt,down_edge);
  }
*/

  int edgCounter = 0;
  for(int edgID=0;edgID<mesh.nElementBoundaries_global;edgID++){
    if(hasModel){
      gEnt = m->findModelEntity(meshBoundary2Model[2*edgID+1],meshBoundary2Model[2*edgID]);
      std::cout<<"What is the search say? "<<m->getModelType(gEnt)<<" "<<m->getModelTag(gEnt)<<" Type "<<meshBoundary2Model[2*edgID+1]<<" ID "<<meshBoundary2Model[2*edgID]<<std::endl;
      if(meshBoundary2Model[2*edgID+1]==1) //if entity is a on a model boundary
        modelBoundaryMaterial[meshBoundary2Model[2*edgID]] = mesh.elementBoundaryMaterialTypes[edgID];
    }
    else{
      if(mesh.exteriorElementBoundariesArray[edgCounter]==edgID){
        gEnt = m->findModelEntity(1,edgCounter);
        modelBoundaryMaterial[edgCounter] = mesh.elementBoundaryMaterialTypes[edgID]; 
        edgCounter++;
      }
      else {
        assert(mesh.elementBoundaryMaterialTypes[edgID]==0);
        gEnt = m->findModelEntity(2,mesh.elementMaterialTypes[0]);
      }
    }
    //modelBoundaryMaterial[edgID] = mesh.elementBoundaryMaterialTypes[edgID]; 
    //int actual_edgID = mesh.exteriorElementBoundariesArray[edgID];
    down_edge[0] = vertices[mesh.edgeNodesArray[edgID*2]];
    down_edge[1] = vertices[mesh.edgeNodesArray[edgID*2+1]];
    //std::cout<<"Mesh Edge Dimension is "<<m->getModelType(gEnt)<<std::endl;
    m->createEntity(apf::Mesh::EDGE,gEnt,down_edge);
  }

  for(int i=0;i<numModelBoundaries;i++){
    std::cout<<"Model boundary "<< i<< " has tag "<< modelBoundaryMaterial[i]<<std::endl;
  }

  std::cout<<"Initializing RECONSTRUCTION!\n";
  modelRegionMaterial[0] = mesh.elementMaterialTypes[0]; //single material type at the moment
  for(int fID=0; fID<mesh.nElements_global;fID++){
      for(int i=0; i<mesh.nNodes_element;i++){
        int vID = mesh.elementNodesArray[fID*mesh.nNodes_element+i];
        elemverts[i]=vertices[vID];
      }
      gEnt = m->findModelEntity(2,mesh.elementMaterialTypes[fID]);
      std::cout<<"Mesh Face Dimension is "<<m->getModelType(gEnt)<<std::endl;
      apf::buildElement(m,gEnt,apf::Mesh::TRIANGLE,elemverts);
  }

  std::cout<<"Looped RECONSTRUCTION!\n";
  m->acceptChanges();
  m->verify();

  gmi_model* g2 = m->getModel();
  gmi_iter* gIter = gmi_begin(g2,0);
  gmi_ent* gmiEnt;
  while(gmiEnt = gmi_next(g2,gIter)){
    std::cout<<"What is the tag? "<<gmi_tag(g2,gmiEnt)<<" material? "<<modelVertexMaterial[gmi_tag(g2,gmiEnt)]<<std::endl;
  }
  gmi_end(g2,gIter);

  //write model and mesh
  gmi_write_dmg(gMod,"Reconstruct.dmg");
  apf::writeVtkFiles("reconstructedMesh", m);

  numberLocally();
  std::cout<<"COMPLETED RECONSTRUCTION!\n";
}

