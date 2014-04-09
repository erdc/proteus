#include "pumi.h"
#include "pumi_mesh.h"
#include "pumi_geom.h"
#include "pumi_geom_geomsim.h"
#include "MeshAdapt.h"
#include "MA.h"
#include "mMeshIO.h"
#include "apfSPR.h"
#include "apfMesh.h"
#include "apfPUMI.h"
#include "maCallback.h"

#include "mpi.h"

#include "MeshAdaptPUMI.h"

MeshAdaptPUMIDrvr::MeshAdaptPUMIDrvr(double Hmax, double Hmin, int NumIter) {
  PUMI_Init(MPI_COMM_WORLD);
  numVar=0;
  hmin=Hmin; hmax=Hmax;
  numIter=NumIter;
  nAdapt=0;
  if(SCUTIL_CommRank()==0)
     printf("Setting hmax=%lf, hmin=%lf, numIters(meshadapt)=%d\n",hmax, hmin, numIter);
}

MeshAdaptPUMIDrvr::~MeshAdaptPUMIDrvr() {
}

int MeshAdaptPUMIDrvr::initProteusMesh(Mesh& mesh) {
  std::cout << "Initializing proteus mesh\n"; 
  // do PUMI stuff to mesh object here.
  return 0;
}

int MeshAdaptPUMIDrvr::readGeomModel(const std::string &geom_sim_file)
{

  PUMI_Geom_RegisterGeomSim(); //commented out for erdc team
  FILE * pFile = fopen(geom_sim_file.c_str(), "r");
  if (pFile)
  {
    std::cerr<<"  saw simmetrix model: "<<geom_sim_file<<"\n";
    PUMI_Geom_LoadFromFile(PUMI_GModel,geom_sim_file.c_str()); //commented out for erdc team
    fclose(pFile);
  }
  else
    PUMI_GModel = 0;

  // print geometry info
  if(PUMI_GModel!=0) { 
    std::cerr<<"[GFace Info] Total Number of Faces of this Model: "<<GM_numFaces(PUMI_GModel)<<std::endl;
    std::cerr<<"[GEdge Info] Total Number of Edges of this Model: "<<GM_numEdges(PUMI_GModel)<<std::endl;
    std::cerr<<"[GVert Info] Total Number of Vertices of this Model: "<<GM_numVertices(PUMI_GModel)<<std::endl;
  }
    
    return 0;
}

int MeshAdaptPUMIDrvr::readPUMIMesh(const char* SMS_fileName){

  PUMI_Mesh_Create(PUMI_GModel, PUMI_MeshInstance);

  if(PCU_Comm_Peers()==1) {
     std::cout << "Reading serial PUMI mesh\n";
     PUMI_Mesh_LoadFromFile(PUMI_MeshInstance,SMS_fileName,0, "pumi"); //0 for serial for now
  } else {
     std::cout << "Reading parallel PUMI mesh\n";
     PUMI_Mesh_LoadFromFile(PUMI_MeshInstance,SMS_fileName,1, "pumi"); //0 for serial for now
  }
  
  int isValid;
  PUMI_Mesh_Verify(PUMI_MeshInstance,&isValid);
  if(isValid) {
    std::cout<<" PUMI mesh verified\n";
  } else {
    std::cout<<" PUMI mesh from file didn't verify, exiting!\n";
    exit(0);
  }
  
  int err = PUMI_Mesh_GetPart(PUMI_MeshInstance, 0, PUMI_Part);

  comm_size = SCUTIL_CommSize();
  comm_rank = SCUTIL_CommRank();
  
  return 0;
}

int MeshAdaptPUMIDrvr::AdaptPUMIMesh() {

  int return_verify;

  PUMI_Mesh_Verify(PUMI_MeshInstance,&return_verify);
  if(return_verify){ std::cerr << "PUMI verify completed succesfully!\n"; } else { exit(0); }
  
  pMAdapt MA_Drvr;
  MA_NewMeshAdaptDrvr_ModelType(MA_Drvr, PUMI_MeshInstance, Application, 2); // third param (0,1,2 : no snapping, non-parametric, parametric)

//  DeleteMeshEntIDs();
  
  apf::Mesh* apf_mesh = apf::createMesh(PUMI_Part);
  getFieldFromTag(apf_mesh, PUMI_MeshInstance,"Solution");
  ma::FieldCallback(MA_Drvr, apf_mesh);
  
  CalculateAnisoSizeField(MA_Drvr, phif);
//  CalculateSizeField(MA_Drvr);
//  CBFunction CB = 0;
  CBFunction CB = TransferTopSCOREC;

  MA_SetCB(MA_Drvr, CB, (void*)this);   // called during mesh modification (see MA.h in meshMeshAdaptPUMI)

  /// Adapt the mesh
  MA_SetNumIt(MA_Drvr, numIter);    // limits the number of iterations of meshMeshAdaptPUMI (splits, collapses, moves...)
  double beta[3] = {2.0, 2.0, 2.0};
  MA_AnisoSmooth(MA_Drvr, beta);
  MA_Adapt(MA_Drvr);  // does the MeshAdaptPUMI
  MA_Del(MA_Drvr);  // deletes the meshMeshAdaptPUMI object
  SCUTIL_Sync();
  
  PUMI_Mesh_DelTag (PUMI_MeshInstance, SFTag, 1);
  if(comm_size>1)
     PUMI_Mesh_DelTag (PUMI_MeshInstance, GlobNumberTag, 1);
  
  PUMI_Mesh_Verify(PUMI_MeshInstance,&return_verify);

//  getTagFromField(apf_mesh, PUMI_MeshInstance, "Solution");

  //partition the mesh
  PUMI_Mesh_SetNumPart(PUMI_MeshInstance, 1);
//  PUMI_Mesh_SetPtnParam(mesh, method, approach, imbTol, dbgLvl);
  PUMI_Mesh_GlobPtn(PUMI_MeshInstance);
  int err = PUMI_Mesh_GetPart(PUMI_MeshInstance, 0, PUMI_Part);

  exportMeshToVTK(PUMI_MeshInstance, "pumi_adapt.vtk");
  if(return_verify){ 
    std::cerr << "Adapted PUMI mesh verify completed succesfully!\n"; 
  } else {
    exit(1);
  }

  nAdapt++; //counter for number of adapt steps

  return 0;
}

