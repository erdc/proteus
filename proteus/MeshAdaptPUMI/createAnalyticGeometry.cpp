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

double boxLength;
double boxWidth;
double boxHeight;
int edgeMap[12] = {50,48,46,52,11,16,20,6,73,72,71,74};
int faceLoop[6] = {80,78,76,82,42,24};


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


void regionFunction(double const p[2], double x[3], void*)
{
  (void)p;
  (void)x;
}


void makeBox(gmi_model* model)
{
  //making a box

  int vertPer = 0;
  double vertRan[1][2] = {{0.0,0.0}};
  int vertexMap[8] = {58,56,54,60,5,10,15,2};
  gmi_ent* g_vert[8];
  g_vert[0] = gmi_add_analytic(model, 0, 58, vert0, &vertPer, vertRan, 0);
  g_vert[1] = gmi_add_analytic(model, 0, 56, vert1, &vertPer, vertRan, 0);
  g_vert[2] = gmi_add_analytic(model, 0, 54, vert2, &vertPer, vertRan, 0);
  g_vert[3] = gmi_add_analytic(model, 0, 60, vert3, &vertPer, vertRan, 0);
  g_vert[4] = gmi_add_analytic(model, 0, 5, vert4, &vertPer, vertRan, 0);
  g_vert[5] = gmi_add_analytic(model, 0, 10, vert5, &vertPer, vertRan, 0);
  g_vert[6] = gmi_add_analytic(model, 0, 15, vert6, &vertPer, vertRan, 0);
  g_vert[7] = gmi_add_analytic(model, 0, 2, vert7, &vertPer, vertRan, 0);


  int edgePer = 0;
  double edgeRan[1][2] = {{0.0,1.0}};
  gmi_ent* g_edge[12];

  g_edge[0] = gmi_add_analytic(model, 1, 50, edge0, &edgePer, edgeRan, 0);
  g_edge[1] = gmi_add_analytic(model, 1, 48, edge1, &edgePer, edgeRan, 0);
  g_edge[2] = gmi_add_analytic(model, 1, 46, edge2, &edgePer, edgeRan, 0);
  g_edge[3] = gmi_add_analytic(model, 1, 52, edge3, &edgePer, edgeRan, 0);
  g_edge[4] = gmi_add_analytic(model, 1, 11, edge4, &edgePer, edgeRan, 0);
  g_edge[5] = gmi_add_analytic(model, 1, 16, edge5, &edgePer, edgeRan, 0);
  g_edge[6] = gmi_add_analytic(model, 1, 20, edge6, &edgePer, edgeRan, 0);
  g_edge[7] = gmi_add_analytic(model, 1, 6, edge7, &edgePer, edgeRan, 0);
  g_edge[8] = gmi_add_analytic(model, 1, 73, edge8, &edgePer, edgeRan, 0);
  g_edge[9] = gmi_add_analytic(model, 1, 72, edge9, &edgePer, edgeRan, 0);
  g_edge[10] = gmi_add_analytic(model, 1, 71, edge10, &edgePer, edgeRan, 0);
  g_edge[11] = gmi_add_analytic(model, 1, 74, edge11, &edgePer, edgeRan, 0);

  //reparameterize vertices on edges
  agm_bdry b;
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

  //make faces

  int facePeriodic[2] = {0, 0};
  double faceRanges[2][2] = {{0,1.0},{0,1.0}};
  gmi_ent* g_face[6];
  g_face[0] = gmi_add_analytic(model, 2, 80, face0, facePeriodic, faceRanges, 0);
  g_face[1] = gmi_add_analytic(model, 2, 78, face1, facePeriodic, faceRanges, 0);
  g_face[2] = gmi_add_analytic(model, 2, 76, face2, facePeriodic, faceRanges, 0);
  g_face[3] = gmi_add_analytic(model, 2, 82, face3, facePeriodic, faceRanges, 0);
  g_face[4] = gmi_add_analytic(model, 2, 42, face4, facePeriodic, faceRanges, 0);
  g_face[5] = gmi_add_analytic(model, 2, 24, face5, facePeriodic, faceRanges, 0);

  //reparam edges onto face

  int edgeLoop[6][4] = {{50,72,11,73},{72,48,71,16},{46,74,20,71},{52,73,6,74},{50,52,46,48},{11,16,20,6}}; 
  int edgeReparamLoop[6][4] = {{0,3,2,1},{1,0,3,2},{0,1,2,3},{0,1,2,3},{0,1,2,3},{0,3,2,1}}; 
  
  typedef void (*ParametricFunctionArray) (double const from[2], double to[2], void*);
  ParametricFunctionArray edgeFaceFunction[] = 
    {
      reparamEdge_0,
      reparamEdge_1,
      reparamEdge_2,
      reparamEdge_3,
    };

  for(int i=0; i<6;i++)
  {
    b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_face[i]));
    for(int j=0; j<4;j++)
    {
      agm_use faceUse = add_adj(model, b, edgeLoop[i][j]);
      gmi_add_analytic_reparam(model, faceUse, edgeFaceFunction[edgeReparamLoop[i][j]], 0);
    }
  }


  gmi_add_analytic_cell(model,3,92);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(gmi_find(model,3,92)));
  for(int i=0; i<6;i++)
  {
    agm_use regionUse = add_adj(model, b, faceLoop[i]);
    gmi_add_analytic_reparam(model, regionUse, regionFunction, 0);
  }

  agm_use regionUse = add_adj(model, b, 123);
  gmi_add_analytic_reparam(model, regionUse, regionFunction, 0);


  return;
}

class Box{
    public:
        void makeBox(gmi_model*);
};

void Box::makeBox(gmi_model* model)
{
  //making a box

  int vertPer = 0;
  double vertRan[1][2] = {{0.0,0.0}};
  int vertexMap[4] = {58,5,10,56};
  gmi_ent* g_vert[4];
  g_vert[0] = gmi_add_analytic(model, 0, 58, vert0, &vertPer, vertRan, 0);
  g_vert[1] = gmi_add_analytic(model, 0, 5, vert1, &vertPer, vertRan, 0);
  g_vert[2] = gmi_add_analytic(model, 0, 10, vert2, &vertPer, vertRan, 0);
  g_vert[3] = gmi_add_analytic(model, 0, 56, vert3, &vertPer, vertRan, 0);

  int edgePer = 0;
  double edgeRan[1][2] = {{0.0,1.0}};
  gmi_ent* g_edge[4];

  g_edge[0] = gmi_add_analytic(model, 1, 50, edge0, &edgePer, edgeRan, 0);
  g_edge[1] = gmi_add_analytic(model, 1, 48, edge1, &edgePer, edgeRan, 0);
  g_edge[2] = gmi_add_analytic(model, 1, 46, edge2, &edgePer, edgeRan, 0);
  g_edge[3] = gmi_add_analytic(model, 1, 52, edge3, &edgePer, edgeRan, 0);

  //reparameterize vertices on edges
  agm_bdry b;
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
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);

  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(g_edge[3]));
  edgeUse0 = add_adj(model, b, vertexMap[3]);
  edgeUse0_1 = add_adj(model,b,vertexMap[0]);
  gmi_add_analytic_reparam(model, edgeUse0, reparamVert_zero, 0);
  gmi_add_analytic_reparam(model, edgeUse0_1, reparamVert_one, 0);

  //reparam edges onto face
  int edgeLoop[4] = {50,48,46,52}; 

  //gmi_add_analytic_cell(model,2,92);
  int faPer[2] = {0, 0};
  double faRan[2][2] = {{0,1},{0,1}};

  typedef void (*ParametricFunctionArray) (double const from[2], double to[2], void*);
  ParametricFunctionArray edgeFaceFunction[] = 
    {
      reparamEdge_0,
      reparamEdge_3,
      reparamEdge_2,
      reparamEdge_1,
    };

  gmi_ent* f = gmi_add_analytic(model, 2, 92, face4, faPer, faRan, 0);
  b = agm_add_bdry(gmi_analytic_topo(model), agm_from_gmi(gmi_find(model,2,92)));
  for(int i=0; i<4;i++)
  {
    agm_use regionUse = add_adj(model, b, edgeLoop[i]);
    gmi_add_analytic_reparam(model, regionUse, edgeFaceFunction[i], 0);
  }

  agm_use regionUse = add_adj(model, b, 123);
  gmi_add_analytic_reparam(model, regionUse, regionFunction, 0);

  return;
}

//create sphere
const double pi = apf::pi;
int sphereFaceID = 123;
double sphereRadius;
double xyz_offset[3];

void sphereFace(double const p[2], double x[3], void*)
{
  x[0] = xyz_offset[0]+sphereRadius*cos(p[0]) * sin(p[1]);
  x[1] = xyz_offset[1]+sphereRadius*sin(p[0]) * sin(p[1]);
  x[2] = xyz_offset[2]+sphereRadius*cos(p[1]);
}

void makeSphere(gmi_model* model)
{
  int faPer[2] = {1, 0};
  double faRan[2][2] = {{0,6.28318530718},{0.0,apf::pi}};

  gmi_add_analytic(model, 2, sphereFaceID, sphereFace, faPer, faRan, 0);
}

class Sphere{
    public:
    int faceID = 123;
    double radius;
    double offset[3];
    int dim; 
    Sphere(int x){
        dim = x;
    }
    //void sphereFaces(double const*, double*, void*);
    //void circleFace(double const, double, void*);
    void makeSphere(gmi_model*);
    
};

void Sphere::makeSphere(gmi_model* model)
{ 
  int faPer[2] = {1, 0};
  double faRan[2][2] = {{0,6.28318530718},{0.0,apf::pi}};
  if(dim==2){
    faRan[1][1] = 0.0; 
  }

  sphereRadius = radius;
  gmi_add_analytic(model, dim-1, sphereFaceID, sphereFace, faPer, faRan, 0);
}

//model tags are based off gmsh default outputs...
void setParameterization(gmi_model* model,apf::Mesh2* m)
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
        m->setModelEntity(ent,(apf::ModelEntity*)gmi_find(model,modelType,sphereFaceID));
      }
      else
        m->setModelEntity(ent,(apf::ModelEntity*)gmi_find(model,modelType,modelTag));

    }
    m->end(it);
  }
  m->setModel(model);
  m->acceptChanges();
    

  //Need to set the parametric coordinates of each of the boundary vertices
  std::map<int,int> edgeParam;
  int edgeScales[12] = {0,1,0,1,0,1,0,1,2,2,2,2};
  double edgeLengths[3] = {boxLength,boxWidth,boxHeight};
  for(int i=0;i<12;i++)
  {
    edgeParam[edgeMap[i]] = edgeScales[i];
  }

  std::map<int,int(*)[2]> faceParam;
  int faceScales[6][2] = {{0,2},{1,2},{0,2},{1,2},{0,1},{0,1}};
  for(int i = 0; i<6;i++)
  {
    faceParam[faceLoop[i]] = &(faceScales[i]); 
  }
  
  apf::MeshIterator* it = m->begin(0);
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
      apf::Vector3 newParam;
      m->getPoint(ent,0,pt);
      m->getParam(ent,oldParam);
      if(modelType==1 && modelTag!=sphereFaceID)
      {
        int relevantIndex = edgeParam[modelTag];
        newParam[0]=pt[relevantIndex]/edgeLengths[relevantIndex];
        m->setParam(ent,newParam);
      }
      else if (modelType==2 && modelTag!=sphereFaceID && m->getDimension()>2)
      {
        int* relevantIndex = faceParam[modelTag][0]; //size is 2
        newParam[0] = pt[relevantIndex[0]]/edgeLengths[relevantIndex[0]];
        newParam[1] = pt[relevantIndex[1]]/edgeLengths[relevantIndex[1]];
        m->setParam(ent,newParam);
      }
      else if (modelType==2 && modelTag == sphereFaceID)
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
      else if (modelType==1 && modelTag == sphereFaceID) //2D version
      {
        double argy = (pt[1]-xyz_offset[1]);
        double argx = (pt[0]-xyz_offset[0]);
        if(argx == 0 && argy ==0)
          newParam[0] = 0.0; // not sure if this will happen or if this is right
        else 
          newParam[0] = atan2(argy,argx);

        m->setParam(ent,newParam);
      }

    } //end if
  } //end while
  m->end(it);
  m->acceptChanges();
}

void setParameterization2D(gmi_model* model,apf::Mesh2* m)
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
      if(modelType==2 && modelTag!=92)
        std::cout<<"model tag,type "<<modelTag<<" "<<modelType<<std::endl;
      if(modelTag > 139 && modelTag< 148)
      {
        m->setModelEntity(ent,(apf::ModelEntity*)gmi_find(model,modelType,sphereFaceID));
      }
      else
        m->setModelEntity(ent,(apf::ModelEntity*)gmi_find(model,modelType,modelTag));

    }
    m->end(it);
  }
  m->setModel(model);
  m->acceptChanges();

  //Need to set the parametric coordinates of each of the boundary vertices
  std::map<int,int> edgeParam;
  const int numEdges = 4;
  int edgeScales[numEdges] = {0,1,0,1};
  double edgeLengths[2] = {boxLength,boxWidth};
  for(int i=0;i<numEdges;i++)
  {
    edgeParam[edgeMap[i]] = edgeScales[i];
  }

  apf::MeshIterator* it = m->begin(0);
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
      if(modelType==1 && modelTag!=sphereFaceID)
      {
        int relevantIndex = edgeParam[modelTag];
        newParam[0]=pt[relevantIndex]/edgeLengths[relevantIndex];
        std::cout<<"model tag "<<modelTag<<" newParam "<<newParam<<" old Param "<<oldParam<<" pt "<<pt<<std::endl;
        m->setParam(ent,newParam);
      }
      else if (modelType==1 && modelTag == sphereFaceID) //2D version
      {
        std::cout<<"xyz offset "<<xyz_offset[0]<<" "<<xyz_offset[1]<<std::endl;
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

void MeshAdaptPUMIDrvr::initialAdapt_analytic(){

  //apf::Field* size_initial = apf::createLagrangeField(m,"size_initial",apf::SCALAR,1);
  size_iso = apf::createLagrangeField(m,"proteus_size",apf::SCALAR,1);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* ent;
  hmin = 0.04;
  hmax = 0.1;
  std::cout<<"hmin, hmax "<<hmin<<" "<<hmax<<std::endl;
  std::cout<<"xyz offset "<<xyz_offset[0]<<" "<<xyz_offset[1]<<" "<<xyz_offset[2]<<std::endl;
  while( (ent = m->iterate(it)) )
  {
/*
    apf::Vector3 pt;
    m->getPoint(ent,0,pt);
    if(sqrt( (pt[0]-xyz_offset[0])*(pt[0]-xyz_offset[0])+ (pt[1]-xyz_offset[1])*(pt[1]-xyz_offset[1]) + (pt[2]-xyz_offset[2])*(pt[2]-xyz_offset[2])) < sphereRadius*1.5)
    {
      apf::setScalar(size_iso,ent,0,hmin);
    }
    else
    {
      apf::setScalar(size_iso,ent,0,hmax);
    }
*/
      apf::setScalar(size_iso,ent,0,hmin);
  }
  m->end(it);

  //gradeMesh(1.5);
   
  std::cout<<"adapt0\n";
  apf::writeVtkFiles("initialProteus",m);
  ma::Input* in = ma::configure(m,size_iso);
  in->maximumIterations = 10;
  in->shouldSnap = true;
  in->shouldTransferParametric = true;
  in->shouldFixShape = true;
  in->debugFolder="./debug_fine";
  ma::adaptVerbose(in,false);
  std::cout<<"adapt1\n";
  m->verify();
  
  //apf::writeVtkFiles("initialAdapt",m);
  freeField(size_iso);

/*
  size_iso = apf::createLagrangeField(m,"proteus_size",apf::SCALAR,1);
  it = m->begin(0);
  while( (ent = m->iterate(it)) )
  {
    apf::Vector3 pt;
    m->getPoint(ent,0,pt);
    if(sqrt( (pt[0]-xyz_offset[0])*(pt[0]-xyz_offset[0])+ (pt[1]-xyz_offset[1])*(pt[1]-xyz_offset[1]) + (pt[2]-xyz_offset[2])*(pt[2]-xyz_offset[2])) < sphereRadius*1.5)
      apf::setScalar(size_iso,ent,0,hmin);
    else
      apf::setScalar(size_iso,ent,0,hmax);
  }
  m->end(it);

  gradeMesh(1.5);
   
  in = ma::configure(m,size_iso);
  in->maximumIterations = 10;
  in->shouldSnap = true;
  in->shouldTransferParametric = true;
  in->shouldFixShape = true;
  in->debugFolder="./debug_fine";
  ma::adaptVerbose(in,false);
  m->verify();
  
  //apf::writeVtkFiles("initialAdapt2",m);
  freeField(size_iso);
*/
}


void MeshAdaptPUMIDrvr::updateSphereCoordinates(double*sphereCenter)
{
  xyz_offset[0] = sphereCenter[0];
  xyz_offset[1] = sphereCenter[1];
  xyz_offset[2] = sphereCenter[2];
}

gmi_model* MeshAdaptPUMIDrvr::createSphereInBox(double* boxDim,double*sphereCenter, double radius)
{
  sphereRadius = radius;
  boxLength = boxDim[0];
  boxWidth = boxDim[1];
  boxHeight = boxDim[2];
  xyz_offset[0] = sphereCenter[0];
  xyz_offset[1] = sphereCenter[1];
  xyz_offset[2] = sphereCenter[2];
  
  //lion_set_verbosity(1);

  //create the analytic model 
  gmi_model* model = gmi_make_analytic();
  
  //add the sphere
 
  makeSphere(model);
  
  //add the box
  makeBox(model);

  //apf::writeVtkFiles("initialInitial",m);
  setParameterization(model,m);
  m->verify();

  //initial adapt
  initialAdapt_analytic();

  m->verify();
  return model;
}

gmi_model* MeshAdaptPUMIDrvr::createCircleInBox(double* boxDim,double*sphereCenter, double radius)
{
  Sphere* circle = new Sphere(2);
  circle->radius = radius;
  boxLength = boxDim[0];
  boxWidth = boxDim[1];
  boxHeight = boxDim[2];
  xyz_offset[0] = sphereCenter[0];
  xyz_offset[1] = sphereCenter[1];
  xyz_offset[2] = sphereCenter[2];
  
  lion_set_verbosity(1);

  //create the analytic model 
  gmi_model* model = gmi_make_analytic();
  
  //add the sphere
 
  //circle->makeSphere(model);
  
  //add the box
  Box* box;  
  box->makeBox(model);

  //apf::writeVtkFiles("initialInitial",m);
  setParameterization2D(model,m); //need to modify
  m->verify();

  //initial adapt
  initialAdapt_analytic();

  m->verify();

  return model;
}

