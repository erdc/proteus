#include "MeshAdaptPUMI.h"
#include <apfMesh.h>

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

    //Rewrite only necessary components
    //Couette 
/*
    m->getPoint(v,0,pt);
    double Lz = 0.05;
    double Uinf = 1.0;
    tmp[0] =0 ; //pressure
    tmp[1] =0; //u
    tmp[2] = Uinf*pt[2]/Lz;
    tmp[3] =0;
*/
    //Poiseuille Flow dpdy=-1
//*
    m->getPoint(v,0,pt);
    double Lz = 0.05;
    double Ly = 0.2;
    tmp[0] = 1-pt[1]/Ly; //pressure starts at 1 and goes to 0
    tmp[1] = 0;
    tmp[2] = 0.5/0.0010021928*(-1)*(pt[2]*pt[2]-Lz*pt[2]); 
    tmp[3] = 0;
//*/

    apf::setComponents(solution, v, 0, &tmp[0]); 
  }
  m->end(it);
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

