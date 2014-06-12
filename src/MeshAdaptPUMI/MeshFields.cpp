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
  assert(nN == static_cast<int>(m->count(0)));
  numVar = nVar;
  solution = apf::createPackedField(m, "proteus_solution", nVar);
  apf::NewArray<double> tmp(nVar);
  apf::MeshEntity* v;
  apf::MeshIterator* it = m->begin(0);
  size_t i = 0;
  while ((v = m->iterate(it))) {
    for(int j = 0; j < nVar; j++)
      tmp[j] = inArray[j * nN + i];
    apf::setComponents(solution, v, 0, &tmp[0]); 
    ++i;
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
  size_t i = 0;
  while ((v = m->iterate(it))) {
    apf::getComponents(solution, v, 0, &tmp[0]); 
    for(int j = 0; j < nVar; j++)
      outArray[j * nN + i] = tmp[j];
    ++i;
  }
  m->end(it);
  apf::destroyField(solution);
}

