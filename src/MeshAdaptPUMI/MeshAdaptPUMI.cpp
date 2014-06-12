#include <gmi_mesh.h>
#include <ma.h>
#include <apfSPR.h>
#include <apfMDS.h>

#include "MeshAdaptPUMI.h"

MeshAdaptPUMIDrvr::MeshAdaptPUMIDrvr(double Hmax, double Hmin, int NumIter) {
  PCU_Comm_Init();
  numVar=0;
  hmin=Hmin; hmax=Hmax;
  numIter=NumIter;
  nAdapt=0;
  if(PCU_Comm_Self()==0)
     printf("Setting hmax=%lf, hmin=%lf, numIters(meshadapt)=%d\n",hmax, hmin, numIter);
  global[0] = global[1] = global[2] = global[3] = 0;
  local[0] = local[1] = local[2] = local[3] = 0;
  solution = 0;
  size_iso = 0;
  size_scale = 0;
  size_frame = 0;
}

MeshAdaptPUMIDrvr::~MeshAdaptPUMIDrvr()
{
  freeField(solution);
  freeField(size_iso);
  freeField(size_scale);
  freeField(size_frame);
}

int loadModelAndMesh(const char* modelFile, const char* meshFile)
{
  m = loadMdsMesh(modelFile, meshFile);
  m->verify();
  comm_size = PCU_Comm_Peers();
  comm_rank = PCU_Comm_Self();
}

struct Anisotropic
{
  Anisotropic(apf::Field* f, apf::Field* s)
  {
    frame = f;
    scale = s;
  }
  apf::Field* frame;
  apf::Field* scale;
  void getValue(Entity* vert, Matrix& r, Vector& h)
  {
    apf::getMatrix(frame, vert, 0, r);
    apf::getVector(scale, vert, 0, h);
  }
};

static ma::Input* configureAnisotropic(apf::Mesh2* m,
    apf::Field* frame, apf::Field* scale)
{
  Anisotropic sf(frame, scale);
  return ma::configure(m, &sf);
}

int MeshAdaptPUMIDrvr::AdaptPUMIMesh()
{
  m->verify();
  
  CalculateAnisoSizeField();

  /// Adapt the mesh
  ma::Input* in = configureAnisotropic(mesh, size_frame, size_scale);
  in->runPreZoltan = true;
  in->runMidDiffusion = true;
  in->runPostZoltan = true;
  in->maximumIterations = numIter;
  ma::adapt(in);
  
  apf::destroyField(sizef);
  if(comm_size>1)
    apf::destroyGlobalNumbering(global);
  
  m->verify();

  apf::writeVtkFiles("pumi_adapt", m);

  nAdapt++; //counter for number of adapt steps

  return 0;
}

