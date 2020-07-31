#ifndef CREATE_ANALYTIC_H
#define CREATE_ANALYTIC_H

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfShape.h>

#include <PCU.h>
#include <gmi_analytic.h>
#include <gmi_mesh.h>
#include <apfMesh2.h>

#include <cassert>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include<map>
#include <math.h>
#include <pcu_util.h>

agm_bdry add_bdry(gmi_model* m, gmi_ent* e);
agm_use add_adj(gmi_model* m, agm_bdry b, int tag);
agm_use add_adj(gmi_model* m, agm_bdry b, int dim, int tag);

extern double boxLength;
extern double boxWidth;
extern double boxHeight;
extern double sphereRadius;
extern double xyz_offset[3];

struct Enclosure{
    
    std::vector<int> vertexMap;
    double vertRan[1][2]={{0.0,0.0}};
    int vertPer=0;

    std::vector<int> edgeMap;
    int edgePer = 0;
    double edgeRan[1][2] = {{0.0,1.0}};

    std::vector<int> faceMap;
    int faPer[2] = {0, 0};
    double faRan[2][2] = {{0,1},{0,1}};

    int regionID = 92; // fixed ID
    void makeBox2D(gmi_model* model); 
    void makeBox3D(gmi_model* model); 
};

class Sphere{
    public:
    int faceID; //values are dictated by spatial tools
    double radius;
    double offset[3];
    int dim; 
    int faPer[2] = {1, 0};
    double faRan[2][2] = {{0,6.28318530718},{0.0,apf::pi}};

    Sphere(int x){
        dim = x;
        if(dim==2){
            faRan[1][1] = 0.0; 
            faceID=6;
        }
        if(dim==3)
            faceID=9;
    }
    Sphere(){}
    void makeSphere(gmi_model* model);
};

namespace Reparam{
    void reparameterizeEntities(gmi_model*model,apf::Mesh2*m,Enclosure box, Sphere sphere);
    void reparameterize2D(gmi_model*model,apf::Mesh2*m,Enclosure box, Sphere sphere);
    void reparameterize3D(gmi_model*model,apf::Mesh2*m,Enclosure box, Sphere sphere);
}

#endif
