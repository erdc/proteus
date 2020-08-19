#include <createAnalyticGeometry.h>
#include "MeshAdaptPUMI.h"
#include <ma.h>
#include <lionPrint.h>

//routines to create analytic sphere in a 3D box

agm_bdry add_bdry(gmi_model* m, gmi_ent* e)
{
  return agm_add_bdry(gmi_analytic_topo(m), agm_from_gmi(e));
}

agm_use add_adj(gmi_model* m, agm_bdry b, int tag)
{
  agm* topo = gmi_analytic_topo(m);
  int dim = agm_dim_from_type(agm_bounds(topo, b).type);
  gmi_ent* de = gmi_find(m, dim - 1, tag);
  return agm_add_use(topo, b, agm_from_gmi(de));
}

//anonymous namespace to make all functions local to this file scope

double xyz_offset[3];
double sphereRadius;
double boxLength;
double boxWidth;
double boxHeight;
int geomDim;
Enclosure modelBox;
Sphere modelSphere;
Sphere modelCircle1;
Sphere modelCircle2;
PiercingCylinder modelPiercingCylinder;

namespace 
{
    void vert0(double const p[2], double x[3], void*)
    {
      (void)p;
      x[0] = 0.0;
      x[1] = 0.0;
      x[2] = 0.0;

    }
    void vert1(double const p[2], double x[3], void*)
    {
      (void)p;
      x[0] = boxLength;
      x[1] = 0.0;
      x[2] = 0.0;
    }

    void vert2(double const p[2], double x[3], void*)
    {
      (void)p;
      x[0] = boxLength;
      x[1] = boxWidth;
      x[2] = 0.0;
    }

    void vert3(double const p[2], double x[3], void*)
    {
      (void)p;
      x[0] = 0.0;
      x[1] = boxWidth;
      x[2] = 0.0;
    }

    void vert4(double const p[2], double x[3], void*)
    {
      (void)p;
      x[0] = 0.0;
      x[1] = 0.0;
      x[2] = boxHeight;

    }
    void vert5(double const p[2], double x[3], void*)
    {
      (void)p;
      x[0] = boxLength;
      x[1] = 0.0;
      x[2] = boxHeight;
    }

    void vert6(double const p[2], double x[3], void*)
    {
      (void)p;
      x[0] = boxLength;
      x[1] = boxWidth;
      x[2] = boxHeight;
    }

    void vert7(double const p[2], double x[3], void*)
    {
      (void)p;
      x[0] = 0.0;
      x[1] = boxWidth;
      x[2] = boxHeight;
    }

    void edge0(double const p[2], double x[3], void*)
    {
      x[0] = boxLength*p[0];
      x[1] = 0.0;
      x[2] = 0.0;
    }

    void edge1(double const p[2], double x[3], void*)
    {
      x[0] = boxLength;
      x[1] = boxWidth*p[0];
      x[2] = 0.0;
    }

    void edge2(double const p[2], double x[3], void*)
    {
      x[0] = boxLength*p[0];
      x[1] = boxWidth;
      x[2] = 0.0;
    }

    void edge3(double const p[2], double x[3], void*)
    {
      x[0] = 0.0;
      x[1] = boxWidth*p[0];
      x[2] = 0.0;
    }

    void edge4(double const p[2], double x[3], void*)
    {
      x[0] = boxLength*p[0];
      x[1] = 0.0;
      x[2] = boxHeight;
    }

    void edge5(double const p[2], double x[3], void*)
    {
      x[0] = boxLength;
      x[1] = boxWidth*p[0];
      x[2] = boxHeight;
    }

    void edge6(double const p[2], double x[3], void*)
    {
      x[0] = boxLength*p[0];
      x[1] = boxWidth;
      x[2] = boxHeight;
    }

    void edge7(double const p[2], double x[3], void*)
    {
      x[0] = 0.0;
      x[1] = boxWidth*p[0];
      x[2] = boxHeight;
    }

    void edge8(double const p[2], double x[3], void*)
    {
      x[0] = 0.0;
      x[1] = 0.0;
      x[2] = boxHeight*p[0];
    }

    void edge9(double const p[2], double x[3], void*)
    {
      x[0] = boxLength;
      x[1] = 0.0;
      x[2] = boxHeight*p[0];
    }

    void edge10(double const p[2], double x[3], void*)
    {
      x[0] = boxLength;
      x[1] = boxWidth;
      x[2] = boxHeight*p[0];
    }

    void edge11(double const p[2], double x[3], void*)
    {
      x[0] = 0.0;
      x[1] = boxWidth;
      x[2] = boxHeight*p[0];
    }

    void face0(double const p[2], double x[3], void*)
    {
      x[0] = boxLength*p[0];
      x[1] = 0.0;
      x[2] = boxHeight*p[1];
    }

    void face1(double const p[2], double x[3], void*)
    {
      x[0] = boxLength;
      x[1] = boxWidth*p[0];
      x[2] = boxHeight*p[1];
    }

    void face2(double const p[2], double x[3], void*)
    {
      x[0] = boxLength*p[0];
      x[1] = boxWidth;
      x[2] = boxHeight*p[1];
    }

    void face3(double const p[2], double x[3], void*)
    {
      x[0] = 0.0;
      x[1] = boxWidth*p[0];
      x[2] = boxHeight*p[1];
    }
    void face4(double const p[2], double x[3], void*)
    {
      x[0] = boxLength*p[0];
      x[1] = boxWidth*p[1];
      x[2] = 0.0;
    }
    void face5(double const p[2], double x[3], void*)
    {
      x[0] = boxLength*p[0];
      x[1] = boxWidth*p[1];
      x[2] = boxHeight;
    }

    void reparamVert_zero(double const from[2], double to[2], void*)
    {
      (void)from;
      to[0] = 0;
      to[1] = 0;
    }
    void reparamVert_one(double const from[2], double to[2], void*)
    {
      (void)from;
      to[0] = 1;
      to[1] = 0;
    }

    //from edge parameterization to face parameterization
    void reparamEdge_0(double const from[2], double to[2], void*)
    {
      to[0] = from[0];
      to[1] = 0.0;
    }

    void reparamEdge_1(double const from[2], double to[2], void*)
    {
      to[0] = 0.0;
      to[1] = 1.0-from[0];
    }
    void reparamEdge_2(double const from[2], double to[2], void*)
    {
      to[0] = 1.0 - from[0];
      to[1] = 1.0;
    }

    void reparamEdge_3(double const from[2], double to[2], void*)
    {
      to[0] = 1.0; 
      to[1] = from[0];
    }

    void reparamREdge_0(double const from[2], double to[2], void*)
    {
      to[0] = from[0];
      to[1] = 0.0;
    }

    void reparamREdge_1(double const from[2], double to[2], void*)
    {
      to[0] = 0.0;
      to[1] = from[0];
    }
    void reparamREdge_2(double const from[2], double to[2], void*)
    {
      to[0] = from[0];
      to[1] = 1.0;
    }

    void reparamREdge_3(double const from[2], double to[2], void*)
    {
      to[0] = 1.0; 
      to[1] = from[0];
    }

    void regionFunction(double const p[2], double x[3], void*)
    {
      (void)p;
      (void)x;
    }

    //from circle to box
    void reparam_Circle(double const from[2], double to[2], void*){

        //given theta, need to get y and x in parameterized form
        double x = sphereRadius*cos(from[0])+xyz_offset[0];
        double y = sphereRadius*sin(from[0])+xyz_offset[1];
        to[0] = x/boxLength;
        to[1] = y/boxWidth;
    }
    void reparam_Circle0(double const from[2], double to[2], void*){

        //given theta, need to get y and x in parameterized form
        double x = sphereRadius*cos(from[0])+xyz_offset[0];
        double y = sphereRadius*sin(from[0])+xyz_offset[1];
        to[0] = x/boxLength;
        to[1] = y/boxWidth;
        to[2] = 0.0;
    }
    void reparam_Circle1(double const from[2], double to[2], void*){

        //given theta, need to get y and x in parameterized form
        double x = sphereRadius*cos(from[0])+xyz_offset[0];
        double y = sphereRadius*sin(from[0])+xyz_offset[1];
        to[0] = x/boxLength;
        to[1] = y/boxWidth;
        to[2] = boxHeight;
    }

    void reparam_CircleCylinder0(double const from[2], double to[2], void*){

        //given theta, need to get y and x in parameterized form
        //double x = sphereRadius*cos(from[0])+xyz_offset[0];
        //double y = sphereRadius*sin(from[0])+xyz_offset[1];
        //to[0] = x/boxLength;
        //to[1] = y/boxWidth;
        to[0] = from[0];
        to[1] = 0.0;//from[1];
        //to[2] = 0;
    }
    void reparam_CircleCylinder1(double const from[2], double to[2], void*){

        //given theta, need to get y and x in parameterized form
        //double x = sphereRadius*cos(from[0])+xyz_offset[0];
        //double y = sphereRadius*sin(from[0])+xyz_offset[1];
        //to[0] = x/boxLength;
        //to[1] = y/boxWidth;
        to[0] = from[0];
        to[1] = 1.0;//from[1];
        //to[2] = 1.0;//boxHeight;
    }


//need to set these functions separately from a member because the function signatures are otherwise modified.
//Likewise the internals need to use static variables to still connect them with other objects

    void sphereFace(double const p[2], double x[3], void*)
    {
      if(geomDim == 2){
        x[0] = xyz_offset[0]+sphereRadius*cos(p[0]);
        x[1] = xyz_offset[1]+sphereRadius*sin(p[0]);
        x[2] = 0.0;
      }
      else if(geomDim == 3){
        x[0] = xyz_offset[0]+sphereRadius*cos(p[0]) * sin(p[1]);
        x[1] = xyz_offset[1]+sphereRadius*sin(p[0]) * sin(p[1]);
        x[2] = xyz_offset[2]+sphereRadius*cos(p[1]);
      }
    }
    void circleFace0(double const p[2], double x[3], void*)
    {
        x[0] = xyz_offset[0]+sphereRadius*cos(p[0]);
        x[1] = xyz_offset[1]+sphereRadius*sin(p[0]);
        x[2] = 0.0;
    }

    void circleFace1(double const p[2], double x[3], void*)
    {
        x[0] = xyz_offset[0]+sphereRadius*cos(p[0]);
        x[1] = xyz_offset[1]+sphereRadius*sin(p[0]);
        x[2] = boxHeight;
    }

    void cylinderFace(double const p[2], double x[3], void*)
    {
        x[0] = xyz_offset[0]+sphereRadius*cos(p[0]);
        //x[1] = xyz_offset[1]+sphereRadius*sin(p[0]) * sin(p[1]);
        x[1] = xyz_offset[1]+sphereRadius*sin(p[0]);
        x[2] = boxHeight*p[1]; //need to double check this
    }


}

void checkEntities(apf::Mesh* m){
  apf::MeshEntity* ent;
  apf::MeshIterator* it;
/*
  for(int i =0;i<4;i++)
  {
    apf::MeshIterator* it = m->begin(i);
    while( (ent = m->iterate(it)))
    {
      apf::ModelEntity* g_ent = m->toModel(ent);
      int modelTag = m->getModelTag(g_ent);
      int modelType = m->getModelType(g_ent);
      std::cout<<"model Tag, type "<<modelTag<<" "<<modelType<<std::endl;
    }
    m->end(it);
  }
*/
  if(m->findField("modelTags"))
    apf::destroyField(m->findField("modelTags"));
  if(m->findField("modelType"))
    apf::destroyField(m->findField("modelType"));

  apf::Field* modelTagge = apf::createLagrangeField(m,"modelTags",apf::SCALAR,1);
  it = m->begin(0);
  while( (ent = m->iterate(it)) )
  {
      apf::ModelEntity* g_ent = m->toModel(ent);
      int modelTag = m->getModelTag(g_ent);
      int modelType = m->getModelType(g_ent);
      apf::setScalar(modelTagge,ent,0,modelTag);
  }
  m->end(it);
  
  apf::Field* modelTyper = apf::createLagrangeField(m,"modelType",apf::SCALAR,1);
  it = m->begin(0);
  while( (ent = m->iterate(it)) )
  {
      apf::ModelEntity* g_ent = m->toModel(ent);
      int modelTag = m->getModelTag(g_ent);
      int modelType = m->getModelType(g_ent);
      apf::setScalar(modelTyper,ent,0,modelType);
  }
  m->end(it);


}

typedef void (*EntityMapArray) (double const p[2], double x[3], void*);
typedef void (*ParametricFunctionArray) (double const from[2], double to[2], void*);


void Enclosure::makeBox2D(gmi_model* model)
{
  //making a box
  vertexMap = {58,5,10,56};
  gmi_ent* g_vert[vertexMap.size()];
  
  EntityMapArray vertexPoints[] = {
        vert0,
        vert1,
        vert2,
        vert3
    };

  for(auto i=0; i<vertexMap.size();i++)
    g_vert[i] = gmi_add_analytic(model, 0, vertexMap[i], vertexPoints[i], &vertPer, vertRan, 0); 

  edgeMap = {1,2,3,4};
  gmi_ent* g_edge[edgeMap.size()];

  EntityMapArray edgeEntities[] = {
        edge0,
        edge1,
        edge2,
        edge3
    };

  for(int i=0;i<edgeMap.size();i++)
    g_edge[i] = gmi_add_analytic(model, 1, edgeMap[i], edgeEntities[i], &edgePer, edgeRan, 0);

  //reparameterize vertices on edges
  agm_bdry b;
  std::vector<std::pair<int,int>> indexPairs = {{0,1}, {1,2}, {2,3}, {3,0}};

  for(int i=0;i<edgeMap.size();i++){

    b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[i]));
    agm_use edgeUse0 = add_adj(model, b, vertexMap[indexPairs[i].first]);
    agm_use edgeUse0_1 = add_adj(model,b,vertexMap[indexPairs[i].second]);
    gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
    gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);
  }

  ParametricFunctionArray edgeFaceFunction[] = 
    {
      reparamEdge_0,
      reparamEdge_3,
      reparamEdge_2,
      reparamEdge_1,
    };

  gmi_ent* f = gmi_add_analytic(model, 2, regionID, face4, faPer, faRan, 0);
  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(gmi_find(model,2,regionID)));
  for(int i=0; i<edgeMap.size();i++)
  {
    agm_use regionUse = add_adj(model, b, edgeMap[i]);
    gmi_add_analytic_reparam(model, regionUse, edgeFaceFunction[i], 0);
  }

  return;
}

void Enclosure::makeBox3D(gmi_model* model)
{
  //making a box

  vertexMap = {58,56,54,60,5,10,15,2};
  gmi_ent* g_vert[vertexMap.size()];
  EntityMapArray vertexPoints[] = {
        vert0,
        vert1,
        vert2,
        vert3,
        vert4,
        vert5,
        vert6,
        vert7
    };

  for(auto i=0; i<vertexMap.size();i++)
    g_vert[i] = gmi_add_analytic(model, 0, vertexMap[i], vertexPoints[i], &vertPer, vertRan, 0); 

  edgeMap = {50,48,46,52,11,16,20,6,73,72,71,74};
  gmi_ent* g_edge[edgeMap.size()];
  EntityMapArray edgeEntities[] = {
        edge0,
        edge1,
        edge2,
        edge3,
        edge4,
        edge5,
        edge6,
        edge7,
        edge8,
        edge9,
        edge10,
        edge11
  };
  for(int i=0;i<edgeMap.size();i++)
    g_edge[i] = gmi_add_analytic(model, 1, edgeMap[i], edgeEntities[i], &edgePer, edgeRan, 0);

  //reparameterize vertices on edges
  agm_bdry b;
  std::vector<std::pair<int,int>> indexPairs = {{0,1}, {1,2}, {2,3}, {3,0},
            {4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}  };


  //std::vector<std::pair<int,int>> indexPairs = {{0,1}, {1,2}, {3,2}, {0,3},
  //          {4,5},{5,6},{7,6},{4,7},{0,4},{1,5},{2,6},{3,7}  };

/*
  for(int i=0;i<edgeMap.size();i++){

    b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[i]));
    agm_use edgeUse0 = add_adj(model, b, vertexMap[indexPairs[i].first]);
    agm_use edgeUse0_1 = add_adj(model,b,vertexMap[indexPairs[i].second]);
    gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
    gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);
  }
*/

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[0]));
  agm_use edgeUse0 = add_adj(model, b, vertexMap[0]);
  agm_use edgeUse0_1 = add_adj(model,b,vertexMap[1]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[1]));
  edgeUse0 = add_adj(model, b, vertexMap[1]);
  edgeUse0_1 = add_adj(model,b,vertexMap[2]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[2]));
  edgeUse0 = add_adj(model, b, vertexMap[2]);
  edgeUse0_1 = add_adj(model,b,vertexMap[3]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_one, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_zero, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[3]));
  edgeUse0 = add_adj(model, b, vertexMap[3]);
  edgeUse0_1 = add_adj(model,b,vertexMap[0]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_one, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_zero, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[4]));
  edgeUse0 = add_adj(model, b, vertexMap[4]);
  edgeUse0_1 = add_adj(model,b,vertexMap[5]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[5]));
  edgeUse0 = add_adj(model, b, vertexMap[5]);
  edgeUse0_1 = add_adj(model,b,vertexMap[6]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[6]));
  edgeUse0 = add_adj(model, b, vertexMap[6]);
  edgeUse0_1 = add_adj(model,b,vertexMap[7]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_one, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_zero, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[7]));
  edgeUse0 = add_adj(model, b, vertexMap[7]);
  edgeUse0_1 = add_adj(model,b,vertexMap[4]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_one, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_zero, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[8]));
  edgeUse0 = add_adj(model, b, vertexMap[0]);
  edgeUse0_1 = add_adj(model,b,vertexMap[4]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[9]));
  edgeUse0 = add_adj(model, b, vertexMap[1]);
  edgeUse0_1 = add_adj(model,b,vertexMap[5]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[10]));
  edgeUse0 = add_adj(model, b, vertexMap[2]);
  edgeUse0_1 = add_adj(model,b,vertexMap[6]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[11]));
  edgeUse0 = add_adj(model, b, vertexMap[3]);
  edgeUse0_1 = add_adj(model,b,vertexMap[7]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);


  //reparameterize vertices on edges
  //make faces

  EntityMapArray faceEntities[] = {
        face0,
        face1,
        face2,
        face3,
        face4,
        face5,
  };

  //faceMap = {80,78,76,82,42,24};
  faceMap = {5,4,3,2,1,6};
  gmi_ent* g_face[faceMap.size()];
  for(int i=0;i<edgeMap.size();i++){
    g_face[i] = gmi_add_analytic(model, 2, faceMap[i], faceEntities[i], faPer, faRan, 0);
  }

  //reparam edges onto face
  int numFaces = faceMap.size();
  int numEdgesFace = 4;

  int edgeLoop[numFaces][numEdgesFace] = {{50,72,11,73},{72,48,71,16},{46,74,20,71},{52,73,6,74},{50,52,46,48},{11,16,20,6}}; 
  int edgeReparamLoop[numFaces][numEdgesFace] = {{0,3,2,1},{1,0,3,2},{0,1,2,3},{0,1,2,3},{0,1,2,3},{0,3,2,1}}; 
  
  //typedef void (*ParametricFunctionArray) (double const from[2], double to[2], void*);
  ParametricFunctionArray edgeFaceFunction[] = 
    {
      reparamREdge_0,
      reparamREdge_1,
      reparamREdge_2,
      reparamREdge_3,
    };

  for(int i=0; i<numFaces;i++)
  {
    b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_face[i]));
    for(int j=0; j<numEdgesFace;j++)
    {
      agm_use faceUse = add_adj(model, b, edgeLoop[i][j]);
      gmi_add_analytic_reparam(model, faceUse, edgeFaceFunction[edgeReparamLoop[i][j]], 0);
    }
  }


  gmi_add_analytic_cell(model,3,regionID);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(gmi_find(model,3,regionID)));
  for(int i=0; i<numFaces;i++)
  {
    agm_use regionUse = add_adj(model, b, faceMap[i]);
    gmi_add_analytic_reparam(model, regionUse, regionFunction, 0);
  }

  return;
}

void PiercingCylinder::makePiercingCylinder(gmi_model* model)
{ 

  sphereRadius = radius;
  gmi_add_analytic(model, dim-1, faceID, cylinderFace, faPer, faRan, 0);
}


void Sphere::makeSphere(gmi_model* model)
{ 

  sphereRadius = radius;
  gmi_add_analytic(model, dim-1, faceID, sphereFace, faPer, faRan, 0);
}


void setParameterizationCylinder(gmi_model* model,apf::Mesh2* m, Enclosure box, PiercingCylinder cylinder,Sphere circle1, Sphere circle2 )
{
  //Get the classification of each entity in the SimMesh
  apf::MeshEntity* ent;
  for(int i =0;i<4;i++)
  {
    apf::MeshIterator* it = m->begin(i);
    while( (ent = m->iterate(it)))
    {
      apf::ModelEntity* g_ent = m->toModel(ent);
      int modelTag = m->getModelTag(g_ent);
      int modelType = m->getModelType(g_ent);
      if(modelTag > 139 && modelTag< 148)
      {
        m->setModelEntity(ent,(apf::ModelEntity*)gmi_find(model,modelType,cylinder.faceID));
      }
      else if(modelTag >= 200){
        //std::cout<<circle1.faceID<<" "<<circle2.faceID<<" "<<modelTag<<std::endl;
        if(modelTag > 200 && modelTag <= 204)
            m->setModelEntity(ent,(apf::ModelEntity*)gmi_find(model,modelType,circle1.faceID));
        else if(modelTag >= 209 && modelTag <= 212)
            m->setModelEntity(ent,(apf::ModelEntity*)gmi_find(model,modelType,circle2.faceID));
      }
      else{
        //std::cout<<"belongs nowhere? "<<modelType<<" "<<modelTag<<std::endl;
        m->setModelEntity(ent,(apf::ModelEntity*)gmi_find(model,modelType,modelTag));
      }

    }
    m->end(it);
  }
  checkEntities(m);
  apf::writeVtkFiles("beforeAccept",m);

  m->setModel(model);
  m->acceptChanges();
  checkEntities(m);
  apf::writeVtkFiles("afterAccept",m);


  Reparam::reparameterizeEntities(model,m,box,cylinder,circle1,circle2);
 
}

//model tags are based off gmsh default outputs...
void setParameterization(gmi_model* model,apf::Mesh2* m, Enclosure box, Sphere sphere)
{

  //Get the classification of each entity in the SimMesh
  apf::MeshEntity* ent;
  for(int i =0;i<4;i++)
  {
    apf::MeshIterator* it = m->begin(i);
    while( (ent = m->iterate(it)))
    {
      apf::ModelEntity* g_ent = m->toModel(ent);
      int modelTag = m->getModelTag(g_ent);
      int modelType = m->getModelType(g_ent);
      if(modelTag > 139 && modelTag< 148)
      {
        m->setModelEntity(ent,(apf::ModelEntity*)gmi_find(model,modelType,sphere.faceID));
      }
      else{
        m->setModelEntity(ent,(apf::ModelEntity*)gmi_find(model,modelType,modelTag));
      }

    }
    m->end(it);
  }
  m->setModel(model);
  m->acceptChanges();
 
  Reparam::reparameterizeEntities(model,m,box,sphere);
 
}

void Reparam::reparameterize3D(gmi_model*model,apf::Mesh2*m,Enclosure box, Sphere sphere){

  //Need to set the parametric coordinates of each of the boundary vertices
  std::map<int,int> edgeParam;
  int edgeScales[12] = {0,1,0,1,0,1,0,1,2,2,2,2};
  double edgeLengths[3] = {boxLength,boxWidth,boxHeight};
  for(int i=0;i<12;i++)
  {
    edgeParam[box.edgeMap[i]] = edgeScales[i];
  }

  std::map<int,int(*)[2]> faceParam;
  int faceScales[6][2] = {{0,2},{1,2},{0,2},{1,2},{0,1},{0,1}};
  for(int i = 0; i<6;i++)
  {
    faceParam[box.faceMap[i]] = &(faceScales[i]); 
  }
  
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* ent;
  while( (ent = m->iterate(it)) )
  {
    apf::ModelEntity* g_ent = m->toModel(ent);
    
    apf::MeshEntity* ev[2];
    m->getDownward(ent,0,ev);
    int modelTag = m->getModelTag(g_ent);
    int modelType = m->getModelType(g_ent);
    if(modelType<3 && modelType!=0)
    {
      apf::Vector3 pt;
      apf::Vector3 newParam;
      m->getPoint(ent,0,pt);
      if(modelType==1)
      {
        int relevantIndex = edgeParam[modelTag];
        newParam[0]=pt[relevantIndex]/edgeLengths[relevantIndex];
        m->setParam(ent,newParam);
        //std::cout<<"modelType? "<<modelType<<" "<<modelTag<<" pt "<<pt<<" "<<newParam<<std::endl;
      }
      else if (modelType==2 && modelTag!=sphere.faceID)
      {
        int* relevantIndex = faceParam[modelTag][0]; //size is 2
        newParam[0] = pt[relevantIndex[0]]/edgeLengths[relevantIndex[0]];
        newParam[1] = pt[relevantIndex[1]]/edgeLengths[relevantIndex[1]];
        m->setParam(ent,newParam);
      }
      else if (modelType==2 && modelTag == sphere.faceID)
      {
        double argy = (pt[1]-xyz_offset[1]);
        double argx = (pt[0]-xyz_offset[0]);
        if(argx == 0 && argy ==0)
          newParam[0] = 0.0; // not sure if this will happen or if this is right
        else 
          newParam[0] = atan2(argy,argx);
        double arg2 = (pt[2]-xyz_offset[2])/sphereRadius;
        if(arg2 < -1.0)
          arg2 = -1.0;
        else if (arg2 > 1.0)
          arg2 = 1.0; 

        newParam[1] = acos(arg2);
        if(newParam[0]<0)
          newParam[0] = newParam[0]+2*apf::pi;
        if(newParam[0]>2*apf::pi)
          newParam[0] = newParam[0]-2*apf::pi;
        
        //this is probably unnecessary
        if(newParam[1]<0.0)
          newParam[1] = -1*newParam[1];
        if(newParam[1]>apf::pi)
          newParam[1] = -1*(newParam[1]-2.0*apf::pi);

        m->setParam(ent,newParam);
      }

    } //end if
  } //end while
  m->end(it);
  m->acceptChanges();

}

void Reparam::reparameterizeEntities(gmi_model*model,apf::Mesh2*m,Enclosure box, Sphere sphere){
    if(m->getDimension()==2)
        reparameterize2D(model,m,box, sphere);
    else
        reparameterize3D(model,m,box, sphere);
}
void Reparam::reparameterizeEntities(gmi_model*model,apf::Mesh2*m,Enclosure box, PiercingCylinder cylinder, Sphere circle1, Sphere circle2){
    //reparameterizeCylinder(model,m,box, cylinder,circle1,circle2);

    //Need to set the parametric coordinates of each of the boundary vertices
    std::map<int,int> edgeParam;
    int edgeScales[12] = {0,1,0,1,0,1,0,1,2,2,2,2};
    double edgeLengths[3] = {boxLength,boxWidth,boxHeight};
    for(int i=0;i<12;i++)
    {
        edgeParam[box.edgeMap[i]] = edgeScales[i];
    }

    std::map<int,int(*)[2]> faceParam;
    int faceScales[6][2] = {{0,2},{1,2},{0,2},{1,2},{0,1},{0,1}};
    for(int i = 0; i<6;i++)
    {
        faceParam[box.faceMap[i]] = &(faceScales[i]); 
    }
  
    apf::MeshIterator* it = m->begin(0);
    apf::MeshEntity* ent;
    while( (ent = m->iterate(it)) )
    {
        apf::ModelEntity* g_ent = m->toModel(ent);
    
        apf::MeshEntity* ev[2];
        m->getDownward(ent,0,ev);
        int modelTag = m->getModelTag(g_ent);
        int modelType = m->getModelType(g_ent);
        if(modelType<3 && modelType!=0)
        {
            apf::Vector3 pt;
            apf::Vector3 newParam;
            m->getPoint(ent,0,pt);
            if(modelType==1)
            {
                if(modelTag == circle1.faceID){
                    double argy = (pt[1]-xyz_offset[1]);
                    double argx = (pt[0]-xyz_offset[0]);
                    if(argx == 0 && argy ==0)
                        newParam[0] = 0.0; // not sure if this will happen or if this is right
                    else 
                        newParam[0] = atan2(argy,argx);
                    m->setParam(ent,newParam);
                }
                else if(modelTag == circle2.faceID){
                    double argy = (pt[1]-xyz_offset[1]);
                    double argx = (pt[0]-xyz_offset[0]);
                    if(argx == 0 && argy ==0)
                        newParam[0] = 0.0; // not sure if this will happen or if this is right
                    else 
                        newParam[0] = atan2(argy,argx);
                    m->setParam(ent,newParam);
                }
                else{
                    int relevantIndex = edgeParam[modelTag];
                    newParam[0]=pt[relevantIndex]/edgeLengths[relevantIndex];
                    m->setParam(ent,newParam);
                }
            }
            else if (modelType==2 && modelTag!=cylinder.faceID)
            { 
                int* relevantIndex = faceParam[modelTag][0]; //size is 2
                newParam[0] = pt[relevantIndex[0]]/edgeLengths[relevantIndex[0]];
                newParam[1] = pt[relevantIndex[1]]/edgeLengths[relevantIndex[1]];
                m->setParam(ent,newParam);
            }
            else if (modelType==2 && modelTag == cylinder.faceID)
            {
                double argy = (pt[1]-xyz_offset[1]);
                double argx = (pt[0]-xyz_offset[0]);
                if(argx == 0 && argy ==0)
                    newParam[0] = 0.0; // not sure if this will happen or if this is right
                else 
                    newParam[0] = atan2(argy,argx);
                newParam[1] = pt[2]/boxHeight;

                m->setParam(ent,newParam);
            }

        } //end if
  } //end while
  m->end(it);
  m->acceptChanges();

}


void Reparam::reparameterize2D(gmi_model*model,apf::Mesh2*m,Enclosure box, Sphere sphere){
  std::map<int,int> edgeParam;
  const int numEdges = 4;
  int edgeScales[numEdges] = {0,1,0,1};
  double edgeLengths[2] = {boxLength,boxWidth};
  for(int i=0;i<numEdges;i++)
  {
    edgeParam[box.edgeMap[i]] = edgeScales[i];
  }

  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* ent;
  while( (ent = m->iterate(it)) )
  {
    apf::ModelEntity* g_ent = m->toModel(ent);
    
    apf::MeshEntity* ev[2];
    m->getDownward(ent,0,ev);
    int modelTag = m->getModelTag(g_ent);
    int modelType = m->getModelType(g_ent);
    if(modelType<3 && modelType!=0)
    {
      apf::Vector3 pt;
      apf::Vector3 oldParam;
      apf::Vector3 newParam(0.0,0.0,0.0);
      m->getPoint(ent,0,pt);
      m->getParam(ent,oldParam);
      if(modelType==1 && modelTag!=sphere.faceID)
      {
        int relevantIndex = edgeParam[modelTag];
        newParam[0]=pt[relevantIndex]/edgeLengths[relevantIndex];
        //std::cout<<"model tag "<<modelTag<<" newParam "<<newParam<<" old Param "<<oldParam<<" pt "<<pt<<std::endl;
        m->setParam(ent,newParam);
      }
      else if (modelType==1 && modelTag == sphere.faceID) //2D version
      {
        //std::cout<<"xyz offset "<<xyz_offset[0]<<" "<<xyz_offset[1]<<std::endl;
        double argy = (pt[1]-xyz_offset[1]);
        double argx = (pt[0]-xyz_offset[0]);
        if(argx == 0 && argy ==0)
          newParam[0] = 0.0; // not sure if this will happen or if this is right
        else 
          newParam[0] = atan2(argy,argx);

        m->setParam(ent,newParam);
      }
      else if (modelType==2)
      {
        newParam[0] = pt[0]/boxLength;
        newParam[1] = pt[1]/boxWidth;
        m->setParam(ent,newParam);
      }

    } //end if
  } //end while
  m->end(it);
  m->acceptChanges();

}

void setParameterization2D(gmi_model* model,apf::Mesh2* m, Enclosure box, Sphere sphere)
{
  //Get the classification of each entity in the SimMesh
  apf::MeshEntity* ent;
  for(int i =0;i<3;i++)
  {
    apf::MeshIterator* it = m->begin(i);
    while( (ent = m->iterate(it)))
    {
      apf::ModelEntity* g_ent = m->toModel(ent);
      int modelTag = m->getModelTag(g_ent);
      int modelType = m->getModelType(g_ent);
      //std::cout<<"modelTag "<<modelTag<<" modelType "<<modelType<<std::endl;
      if(modelTag > 139 && modelTag< 148)
      {
        m->setModelEntity(ent,(apf::ModelEntity*)gmi_find(model,modelType,sphere.faceID));
      }
      else
        m->setModelEntity(ent,(apf::ModelEntity*)gmi_find(model,modelType,modelTag));

    }
    m->end(it);
  }
  m->setModel(model);
  m->acceptChanges();

  //Need to set the parametric coordinates of each of the boundary vertices
  Reparam::reparameterizeEntities(model,m,box,sphere);
}

void MeshAdaptPUMIDrvr::initialAdapt_analytic(){
  //at this point, hmin and hmax haven't been set yet
  //apf::Field* size_initial = apf::createLagrangeField(m,"size_initial",apf::SCALAR,1);
  lion_set_verbosity(1);
  apf::Field* size_init = apf::createLagrangeField(m,"proteus_init",apf::SCALAR,1);
  sizeFieldList.push(size_init);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* ent;
  hmin = 0.0125;
  hmax = 0.2;
  while( (ent = m->iterate(it)) )
  {
    apf::Vector3 pt;
    m->getPoint(ent,0,pt);
/*
    if(sqrt( (pt[0]-xyz_offset[0])*(pt[0]-xyz_offset[0])+ (pt[1]-xyz_offset[1])*(pt[1]-xyz_offset[1]) + (pt[2]-xyz_offset[2])*(pt[2]-xyz_offset[2])) < sphereRadius*1.5)
    {
      apf::setScalar(size_iso,ent,0,hmin);
    }
    else
    {
      apf::setScalar(size_iso,ent,0,hmax);
    }
*/
      //apf::setScalar(size_iso,ent,0,0.2);
      apf::ModelEntity* g_ent = m->toModel(ent);
      int modelTag = m->getModelTag(g_ent);
      int modelType = m->getModelType(g_ent);
      if(modelTag == modelPiercingCylinder.faceID && modelType==2 || modelType==1 && (modelTag == modelCircle1.faceID || modelTag == modelCircle2.faceID))
          apf::setScalar(size_init,ent,0,0.0125);
      else
          apf::setScalar(size_init,ent,0,0.1025);
  }
  m->end(it);

  apf::writeVtkFiles("pregradeProteus",m);
  std::cout<<"grade mesh intiial\n";
  //gradeMesh(1.5);
  isotropicIntersect();
  std::cout<<"grade mesh post\n";
  
  checkEntities(m);
  apf::writeVtkFiles("initialProteus",m);
  ma::Input* in = ma::configure(m,size_init);
  in->maximumIterations = 10;
  in->shouldSnap = true;
  in->shouldTransferParametric = true;
  in->shouldFixShape = true;
  in->debugFolder="./debug_fine";
  ma::adaptVerbose(in,true);
  m->verify();
  apf::writeVtkFiles("middleProteus",m);
  
  //apf::writeVtkFiles("initialAdapt",m);
  freeField(size_init);

  size_init = apf::createLagrangeField(m,"proteus_initial",apf::SCALAR,1);
  sizeFieldList.push(size_init);
  it = m->begin(0);
  while( (ent = m->iterate(it)) )
  {
/*
    apf::Vector3 pt;
    m->getPoint(ent,0,pt);
    if(sqrt( (pt[0]-xyz_offset[0])*(pt[0]-xyz_offset[0])+ (pt[1]-xyz_offset[1])*(pt[1]-xyz_offset[1]) + (pt[2]-xyz_offset[2])*(pt[2]-xyz_offset[2])) < sphereRadius*1.5)
      apf::setScalar(size_iso,ent,0,hmin);
    else
      apf::setScalar(size_iso,ent,0,hmax);
*/
      //apf::setScalar(size_iso,ent,0,0.05);
      apf::ModelEntity* g_ent = m->toModel(ent);
      int modelTag = m->getModelTag(g_ent);
      int modelType = m->getModelType(g_ent);
      if(modelTag == modelPiercingCylinder.faceID && modelType==2 || modelType==1 && (modelTag == modelCircle1.faceID || modelTag == modelCircle2.faceID))
          apf::setScalar(size_init,ent,0,0.0125);
      else
          apf::setScalar(size_init,ent,0,0.1025);
    
  }
  m->end(it);

  //gradeMesh(1.5);
   
  in = ma::configure(m,size_init);
  in->maximumIterations = 10;
  in->shouldSnap = true;
  in->shouldTransferParametric = true;
  in->shouldFixShape = true;
  in->debugFolder="./debug_fine2";
  ma::adaptVerbose(in,true);
  m->verify();
  
  apf::writeVtkFiles("finalProteus",m);
  //apf::writeVtkFiles("initialAdapt2",m);
  freeField(size_init);
}


void MeshAdaptPUMIDrvr::updateSphereCoordinates(double*sphereCenter)
{
  xyz_offset[0] = sphereCenter[0];
  xyz_offset[1] = sphereCenter[1];
  xyz_offset[2] = sphereCenter[2];

  char buffer[100];
  sprintf(buffer,"Checking coordinates at update %f %f %f",xyz_offset[0],xyz_offset[1],xyz_offset[2]);
  logEvent(buffer,3);

}

int Splitter::partitionFactor;
void MeshAdaptPUMIDrvr::createAnalyticGeometry(int dim, double* boxDim,double*sphereCenter, double radius)
{
  boxLength = boxDim[0];
  boxWidth = boxDim[1];
  boxHeight = boxDim[2];
  updateSphereCoordinates(sphereCenter);

  //create analytic model
  gmi_model* model = gmi_make_analytic();
  
  //add sphere

  Sphere sphere = Sphere(dim);
  sphere.radius = radius;

  sphere.makeSphere(model);
  
  //add the box
  geomDim = dim;
  Enclosure box;
  if(dim==3){
      box.makeBox3D(model);
      agm_bdry b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(gmi_find(model,dim,box.regionID)));
      agm_use regionUse = add_adj(model, b, sphere.faceID);
      gmi_add_analytic_reparam(model, regionUse, regionFunction, 0);
      setParameterization(model,m,box,sphere);
  }
  else{
      box.makeBox2D(model);  
      agm_bdry b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(gmi_find(model,dim,box.regionID)));
      agm_use regionUse = add_adj(model, b, sphere.faceID);
      gmi_add_analytic_reparam(model, regionUse, reparam_Circle, 0);
      setParameterization2D(model,m,box,sphere);
  }

  modelBox = box;
  modelSphere = sphere;
  isAnalytic = 1;
  m->verify();

  //m->writeNative("analyticMesh.smb");
  apf::Migration* plan = 0;
  gmi_register_mesh();
  Splitter::partitionFactor = PCU_Comm_Peers();
  Splitter::switchToOriginals();
  bool isOriginal = ((PCU_Comm_Self() % Splitter::partitionFactor) == 0);
  if (isOriginal) {
    plan = Splitter::getPlan(m);
  }
  Splitter::switchToAll();
  m = apf::repeatMdsMesh(m, model, plan, Splitter::partitionFactor);
  Parma_PrintPtnStats(m, "");
  //m->writeNative("analyticMesh.smb");

  //initialAdapt_analytic();
  //initialAdapt_analytic(1.0);
  //initialAdapt_analytic(0.1);

  return;

}


void MeshAdaptPUMIDrvr::createAnalyticGeometryCylinder(int dim, double* boxDim,double*sphereCenter, double radius)
{
  boxLength = boxDim[0];
  boxWidth = boxDim[1];
  boxHeight = boxDim[2];
  updateSphereCoordinates(sphereCenter);

  //create analytic model
  gmi_model* model = gmi_make_analytic();
  
  //add sphere

  PiercingCylinder cylinder = PiercingCylinder();
  cylinder.radius = radius;

  cylinder.makePiercingCylinder(model);
  //int edgeID[2] = {200,201};
  int edgeID[2] = {9,10};
  Sphere circle1 = Sphere(2);
  circle1.faceID = edgeID[0];
  circle1.radius = radius;
  //can't use makeSphere because of paramterization function is not general
  gmi_add_analytic(model, 1, circle1.faceID, circleFace0, circle1.faPer, circle1.faRan, 0);

  Sphere circle2 = Sphere(2);
  circle2.faceID = edgeID[1];
  circle2.radius = radius;
  //circle2.makeSphere(model);
  gmi_add_analytic(model, 1, circle2.faceID, circleFace1, circle2.faPer, circle2.faRan, 0);

  //add the box
  geomDim = dim;
  Enclosure box;
  
  //add special flag to makeBox to include circles? 
  //need to classify vertices on these circles independently from the cylindrical face
  
  box.makeBox3D(model);
  //create 2 circles, parameterize it as part of the two faces 1 and 6, faces 4+5

  //face ids come from faceMap
  agm_bdry b;
  agm_use regionUse;

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(gmi_find(model,2,1)));
  regionUse = add_adj(model, b, circle1.faceID);
  gmi_add_analytic_reparam(model, regionUse, reparam_Circle0, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(gmi_find(model,2,6)));
  regionUse = add_adj(model, b, circle2.faceID);
  gmi_add_analytic_reparam(model, regionUse, reparam_Circle1, 0);

  //parameterize the cylinder to the region
  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(gmi_find(model,dim,box.regionID)));
  regionUse = add_adj(model, b, cylinder.faceID);
  gmi_add_analytic_reparam(model, regionUse, regionFunction, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(gmi_find(model,2,cylinder.faceID)));
  regionUse = add_adj(model, b, circle1.faceID);
  gmi_add_analytic_reparam(model, regionUse, reparam_CircleCylinder0, 0);
  regionUse = add_adj(model, b, circle2.faceID);
  gmi_add_analytic_reparam(model, regionUse, reparam_CircleCylinder1, 0);

  setParameterizationCylinder(model,m,box,cylinder,circle1,circle2);

  modelBox = box;
  modelCircle1 = circle1;
  modelCircle2 = circle2;
  modelPiercingCylinder = cylinder;
  isAnalytic = 1;
  m->verify();

  apf::Migration* plan = 0;
  gmi_register_mesh();
  Splitter::partitionFactor = PCU_Comm_Peers();
  Splitter::switchToOriginals();
  bool isOriginal = ((PCU_Comm_Self() % Splitter::partitionFactor) == 0);
  if (isOriginal) {
    plan = Splitter::getPlan(m);
  }
  Splitter::switchToAll();
  m = apf::repeatMdsMesh(m, model, plan, Splitter::partitionFactor);
  Parma_PrintPtnStats(m, "");

  return;

}


void Splitter::freeMesh(apf::Mesh* m)
{
    m->destroyNative();
    apf::destroyMesh(m);
}

apf::Migration* Splitter::getPlan(apf::Mesh* m)
{
/*
    apf::Splitter* splitter = apf::makeZoltanSplitter(
        m, apf::GRAPH, apf::PARTITION, false);
    apf::MeshTag* weights = Parma_WeighByMemory(m);
    apf::Migration* plan = splitter->split(weights, 1.05, partitionFactor);
*/

    apf::Splitter* splitter = Parma_MakeRibSplitter(m);
    apf::MeshTag* weights = Parma_WeighByMemory(m); 
    apf::Migration* plan = splitter->split(weights, 1.10, partitionFactor);


    apf::removeTagFromDimension(m, weights, m->getDimension());
    m->destroyTag(weights);
    delete splitter;
    return plan;
}

void Splitter::switchToOriginals()
{
    int self = PCU_Comm_Self();
    int groupRank = self / Splitter::partitionFactor;
    int group = self % Splitter::partitionFactor;
    MPI_Comm groupComm;
    MPI_Comm_split(MPI_COMM_WORLD, group, groupRank, &groupComm);
    PCU_Switch_Comm(groupComm);
}

void Splitter::switchToAll()
{
    MPI_Comm prevComm = PCU_Get_Comm();
    PCU_Switch_Comm(MPI_COMM_WORLD);
    MPI_Comm_free(&prevComm);
    PCU_Barrier();
}


