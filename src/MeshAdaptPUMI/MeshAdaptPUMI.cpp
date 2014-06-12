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
}

MeshAdaptPUMIDrvr::~MeshAdaptPUMIDrvr() {
}

int MeshAdaptPUMIDrvr::initProteusMesh(Mesh& mesh) {
  /* consider deleting this function */
  return 0;
}

int loadModelAndMesh(const char* modelFile, const char* meshFile)
{
  m = loadMdsMesh(modelFile, meshFile);
  m->verify();
  comm_size = PCU_Comm_Peers();
  comm_rank = PCU_Comm_Self();
}

int MeshAdaptPUMIDrvr::AdaptPUMIMesh()
{
  m->verify();
  
  CalculateAnisoSizeField(MA_Drvr, phif);

  /// Adapt the mesh
  ma::Input* in = ma::configure(m, sizef);
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

