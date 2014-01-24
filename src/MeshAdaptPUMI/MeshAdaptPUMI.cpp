#include "pumi.h"
#include "pumi_mesh.h"
//#include "MeshAdapt.h"
//#include "AcisModel.h"

#include "mpi.h"

#include "MeshAdaptPUMI.h"

//int PUMI_Geom_CreateFromAcisFile(pGModel & model,const char * fileName,const char * tagName);

MeshAdaptPUMIDrvr::MeshAdaptPUMIDrvr() {
}

MeshAdaptPUMIDrvr::~MeshAdaptPUMIDrvr() {
}

int MeshAdaptPUMIDrvr::initProteusMesh(Mesh& mesh) {
  std::cout << "Initializing proteus mesh\n"; 
  initializeMesh(mesh);
  return 0;
}

int MeshAdaptPUMIDrvr::deleteProteusMesh(Mesh& mesh) {
  std::cout << "Deleting proteus mesh\n";
  deleteMesh(mesh);
}

int MeshAdaptPUMIDrvr::readGeomModel(const std::string &acis_geom_file_name)
{
  PUMI_Init(MPI_COMM_WORLD);
  FILE * pFile = fopen(acis_geom_file_name.c_str(), "r");
  if (pFile)
  {
    std::cerr<<"  saw acis model: "<<acis_geom_file_name<<"\n";
    fclose(pFile);
    const char* SCOREC_TAGNAME = "SCORECTAG";
//    PUMI_Geom_CreateFromAcisFile(PUMI_GModel,acis_geom_file_name.c_str(),SCOREC_TAGNAME); //create ACIS model from file
  }
  else
    PUMI_GModel = 0;

  // print geometry info 
//    std::cerr<<"[GFace Info] Total Number of Faces of this Model: "<<GM_numFaces(PUMI_GModel)<<std::endl;
//    std::cerr<<"[GEdge Info] Total Number of Edges of this Model: "<<GM_numEdges(PUMI_GModel)<<std::endl;
//    std::cerr<<"[GVert Info] Total Number of Vertices of this Model: "<<GM_numVertices(PUMI_GModel)<<std::endl;
    
    return 0;
}

int MeshAdaptPUMIDrvr::readPUMIMesh(const char* SMS_fileName){

  PUMI_Mesh_Create(PUMI_GModel, PUMI_MeshInstance);

  if(PCU_Comm_Peers()==1)
     PUMI_Mesh_LoadFromFile(PUMI_MeshInstance,SMS_fileName,0, "pumi"); //0 for serial for now
  else
     PUMI_Mesh_LoadFromFile(PUMI_MeshInstance,SMS_fileName,1, "pumi"); //0 for serial for now
/*
  std::vector<pPart> meshes;
  int rank;
  FMDB_Mesh_GetAllPart(PUMI_MeshInstance, SCUTIL_CommRank(), meshes);
  pMesh mesh = meshes.at(0);

  Mesh_InitBLs(mesh, PUMI_GModel);
*/

  int isValid;
  PUMI_Mesh_Verify(PUMI_MeshInstance,&isValid);
  if(isValid) {
    std::cout<<" PUMI mesh verified\n";
  } else {
    std::cout<<" PUMI mesh from file didn't verify, exiting!\n";
    exit(0);
  }
  std::cout<<" loading PUMI mesh from files\n";

  PUMI_Mesh_GetPart(PUMI_MeshInstance,0,PUMI_Part);
  int err = PUMI_Mesh_GetPart(PUMI_MeshInstance, 0, PUMI_Part);
  
  return 0;
}

int MeshAdaptPUMIDrvr::MeshAdaptPUMI() {

  int return_verify;

  PUMI_Mesh_Verify(PUMI_MeshInstance,&return_verify);
  if(return_verify){ std::cerr << "FMDB verify completed succesfully!\n"; } else { exit(0); }

//  MA_NewMeshMeshAdaptPUMIDrvr_ModelType(MA_Drvr, PUMI_MeshInstance, Application, 2); // third param (0,1,2 : no snapping, non-parametric, parametric)
 // set up size field on verts below
//  PUMI_Mesh_FindTag(PUMI_MeshInstance, "NodeMeshSize", PUMI_SFTag); //todo: replace tags with apf
  double dVtxSize;
  pPartEntIter EntIt;
  pMeshEnt meshEnt;
  PUMI_PartEntIter_Init (PUMI_Part, PUMI_VERTEX, PUMI_POINT, EntIt);
  int isEnd = 0;
  int sizeCounter = 0;

  while (!isEnd)
  {
    PUMI_PartEntIter_GetNext(EntIt, meshEnt);
//    if(SCUtil_SUCCESS == PUMI_MeshEnt_GetDblTag (PUMI_MeshInstance, meshEnt, PUMI_SFTag, &dVtxSize)){ todo: replace with apf
//      MA_SetIsoVtxSize(MA_Drvr, (pVertex)meshEnt, dVtxSize);   // sets size field from tag data
      sizeCounter++;
//    }
    PUMI_PartEntIter_IsEnd(EntIt, &isEnd);
  }
  std::cerr<<" - set size field for "<<sizeCounter<<" vertices\n";

  PUMI_PartEntIter_Del(EntIt);
  isEnd = 0;

//  CBFunction CB = 0;
//  CBFunction CB = TransferTopSCOREC;

//  MA_SetCB(MA_Drvr, CB, (void*)this);   // called during mesh modification (see MA.h in meshMeshAdaptPUMI)

  /// MeshAdaptPUMI the mesh
//  MA_SetNumIt(MA_Drvr, 3);    // limits the number of iterations of meshMeshAdaptPUMI (splits, collapses, moves...)

//  MA_MeshAdaptPUMI(MA_Drvr);  // does the MeshAdaptPUMI

//  MA_Del(MA_Drvr);  // deletes the meshMeshAdaptPUMI object

  SCUTIL_Sync();

  PUMI_Mesh_Verify(PUMI_MeshInstance,&return_verify);
  if(return_verify){ std::cerr << "FMDB verify completed succesfully!\n"; }

  return 0;
}

