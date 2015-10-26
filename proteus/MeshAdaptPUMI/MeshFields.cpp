#include "MeshAdaptPUMI.h"
#include <apfMesh.h>
#include <apfShape.h>

void MeshAdaptPUMIDrvr::freeField(apf::Field*& f)
{
  if (f) {
    apf::destroyField(f);
    f = 0;
  }
}

void MeshAdaptPUMIDrvr::freeNumbering(apf::Numbering*& n)
{
  if (n) {
    apf::destroyNumbering(n);
    n = 0;
  }
}

int MeshAdaptPUMIDrvr::TransferSolutionToPUMI(double* inArray, int nVar, int nN)
{
  assert(nN == static_cast<int>(m->count(0)));
  numVar = nVar;
  solution = apf::createPackedField(m, "proteus_solution", nVar);
  apf::NewArray<double> tmp(nVar);
  apf::MeshEntity* v;
  apf::MeshIterator* it = m->begin(0);
  apf::Vector3 pt;
//std::cout<<"Solution being overwritten "<<std::endl;
  while ((v = m->iterate(it))) {
    int i = localNumber(v);
    for(int j = 0; j < nVar; j++)
      tmp[j] = inArray[j * nN + i];

    int casenum = 0;
    //Rewrite only necessary components
    if(casenum ==0){
    //Poiseuille Flow dpdy=-1
            m->getPoint(v,0,pt);
            double Lz = 0.05;
            double Ly = 0.2;
            tmp[0] = 1-pt[1]/Ly; //pressure starts at 1 and goes to 0
            tmp[1] = 0;
            tmp[2] = 0.5/0.0010021928*(-1/Ly)*(pt[2]*pt[2]-Lz*pt[2]);  //dpdy = 1/Ly
            tmp[3] = 0;
    }
    else if(casenum==1){
    //Couette 
            m->getPoint(v,0,pt);
            double Lz = 0.05;
            double Uinf = 1.0;//2e-3;
            tmp[0] =0 ; //pressure
            tmp[1] =0; //u
            tmp[2] = Uinf*pt[2]/Lz;
            tmp[3] =0;
    }
    apf::setComponents(solution, v, 0, &tmp[0]); 
  }
  m->end(it);
  return 0;
}

int MeshAdaptPUMIDrvr::TransferPropertiesToPUMI(double* rho_p, double* nu_p)
{ 
 rho[0] = rho_p[0]; rho[1] = rho_p[1];
 nu[0] = nu_p[0]; nu[1] = nu_p[1];
 return 0;
}



int MeshAdaptPUMIDrvr::TransferSolutionToProteus(double* outArray, int nVar, int nN)
{
  assert(nN == static_cast<int>(m->count(0)));
  apf::NewArray<double> tmp(nVar);
  apf::MeshEntity* v;
  apf::MeshIterator* it = m->begin(0);
  while ((v = m->iterate(it))) {
    int i = localNumber(v);
    apf::getComponents(solution, v, 0, &tmp[0]); 
    for(int j = 0; j < nVar; j++)
      outArray[j * nN + i] = tmp[j];
  }
  m->end(it);
  apf::destroyField(solution);
  return 0;
}

int MeshAdaptPUMIDrvr::TransferBCtagsToProteus()
{
  //Suppose I have a list of identifiers from Proteus that classifies each boundary element
  apf::MeshIterator* it= m->begin(0);
  apf::MeshEntity* f;
  apf::ModelEntity* me;
  apf::ModelEntity* boundary_face; 
  int tag = 0;
  BCtag[0] =  m->createIntTag("BCtype_p",1);
  BCtag[1] =  m->createIntTag("BCtype_u",1);
  BCtag[2] =  m->createIntTag("BCtype_v",1);
  BCtag[3] =  m->createIntTag("BCtype_w",1);
  while(f=m->iterate(it)){
    int i = localNumber(f);
    me=m->toModel(f);
    tag = m->getModelTag(me);
    boundary_face = m->findModelEntity(2,tag); //faces
    if(me==boundary_face){ //is on model entity
      //Assign a tag to the face for the given type of boundary condition
      int type = 1; //type None:0 DBC:1, NBC:2
      m->setIntTag(f,BCtag[0],&type);
      m->setIntTag(f,BCtag[1],&type);
      m->setIntTag(f,BCtag[2],&type);
      m->setIntTag(f,BCtag[3],&type);
    }
  }
  m->end(it);
  std::cout<<"Finished Transfer of BC Tags "<<std::endl;
  return 0;
}

int MeshAdaptPUMIDrvr::TransferBCsToProteus()
{
  //Want to use some sort of Hierarchic projection 
  apf::FieldShape* BC_shape = apf::getHierarchic(2);
  apf::MeshEntity* v;  
  int nVar = 4; //pressure + 3 velocity components
  //DBC = apf::createPackedField(m, "proteus_DBC", nVar);
  fluxBC = apf::createField(m, "proteus_fluxBC",apf::VECTOR, BC_shape);

  //iterate through faces, find adjacent vertices and edges and then construct the linear system to solve for the coefficients
  apf::MeshIterator* it= m->begin(2);
  apf::MeshEntity* f;
  apf::ModelEntity* me;
  apf::ModelEntity* boundary_face; 
  int tag = 0;
  int count = 0;
  int numqpt;
  apf::Adjacent adjvtx, adjedg;
  while(f=m->iterate(it)){

    if(count==0){ //happens only once
      apf::MeshElement* sample_elem = apf::createMeshElement(m,f);
      numqpt = apf::countIntPoints(sample_elem,integration_order);
      fluxtag[1]= m->createDoubleTag("u_flux",numqpt);
      fluxtag[2]= m->createDoubleTag("v_flux",numqpt);
      fluxtag[3]= m->createDoubleTag("w_flux",numqpt);
      apf::destroyMeshElement(sample_elem);
    }

    me=m->toModel(f);
    tag=m->getModelTag(me);
    boundary_face = m->findModelEntity(2,tag); //faces
    if(me==boundary_face){
      double test[numqpt];
      for(int i=0;i<numqpt;i++){test[i]=i;} 
      m->setDoubleTag(f,fluxtag[1],test);
      double data[numqpt];
      m->getDoubleTag(f,fluxtag[1],&data[0]);
//std::cout<<"What is the data? "<<data[0]<<" "<<data[5]<<" numqpt "<<numqpt<<std::endl;
      //get adjacent vertices
      m->getAdjacent(f,0,adjvtx); 
      //evaluate boundary condition and set value

      //get adjacent edges 
      m->getAdjacent(f,1,adjedg);
    }
  }
  std::cout<<"Finished Transferring BC "<<std::endl;
  m->end(it);
  return 0;
}
