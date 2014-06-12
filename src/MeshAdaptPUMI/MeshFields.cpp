#include "MeshAdaptPUMI.h"
#include <apfSPR.h>
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
  assert(nN == mesh->count(0));

  numVar = nVar;
  solution = apf::createPackedField(mesh, "proteus_solution", nVar);

  apf::NewArray<double> tmp(nVar);
  apf::MeshEntity* v;
  apf::MeshIterator* it = m->begin(0);
  size_t i = 0;
  while ((v = m->iterate(it))) {
    for(int j = 0; j < nVar; j++)
      tmp[j] = inArray[j * nN + i];
    apf::setComponents(solution, v, 0, tmp); 
    ++i;
  }
  m->end(it);
}

int MeshAdaptPUMIDrvr::TransferSolutionToProteus(double* outArray, int nVar, int nN)
{
  assert(nN == mesh->count(0));
  apf::NewArray<double> tmp(nVar);
  apf::MeshEntity* v;
  apf::MeshIterator* it = m->begin(0);
  size_t i = 0;
  while ((v = m->iterate(it))) {
    apf::getComponents(solution, v, 0, tmp); 
    for(int j = 0; j < nVar; j++)
      inArray[j * nN + i] = tmp[j];
    ++i;
  }
  m->end(it);
  apf::destroyField(solution);
}

