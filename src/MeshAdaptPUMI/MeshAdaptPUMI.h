#include <cstdlib>
//#include "MA.h"
#include "pumi_mesh.h"
#include "pumi.h"
#include "pumi_geom.h"
#include "mesh.h"
#include "cmeshToolsModule.h"

class MeshAdaptPUMIDrvr{
 
  public:
  MeshAdaptPUMIDrvr(); 
  ~MeshAdaptPUMIDrvr();

  Mesh mesh_proteus;
  int initProteusMesh(Mesh& mesh);

  int readGeomModel(const std::string &acis_geom_file_name);
  int readPUMIMesh(const char* SMS_fileName);
  int helloworld(const char* hello) { std::cout << hello << "\n"; return 0;}

  //Functions to construct proteus mesh data structures
  int ConstructFromSerialPUMIMesh(Mesh& mesh);
  int ConstructFromParallelPUMIMesh(Mesh& mesh, Mesh& subdomain_mesh);

  int UpdateMaterialArrays(Mesh& mesh, int bdryID, int GeomTag);
  int MeshAdaptPUMI();

  private: 
  pMeshMdl PUMI_MeshInstance;
  pGeomMdl PUMI_GModel;
  std::vector<pPart> PUMI_Parts;
  pPart PUMI_Part;
  int comm_size, comm_rank;
  int elms_owned, faces_owned, edges_owned, vtx_owned;

  pTag elementGlobNumberTag, nodeGlobNumberTag, faceGlobNumberTag, edgeGlobNumberTag;
  pTag GlobNumberTag;
//  pMMeshAdaptPUMI MA_Drvr;

  int ConstructGlobalNumbering(Mesh& mesh);
  int ConstructGlobalStructures(Mesh& mesh);

  int ConstructElements(Mesh& mesh);
  int ConstructNodes(Mesh& mesh);
  int ConstructBoundaries(Mesh& mesh);
  int ConstructEdges(Mesh& mesh);
  int ConstructMaterialArrays(Mesh& mesh);
  
  int CalculateOwnedEnts(PUMI_EntType EntType, int &nOwnedEnt);
  int CommunicateOwnedNumbers(int toSend, int *toReceive);
  int SetOwnerGlobNumbering(pTag, PUMI_EntType, int);
  int SetCopyGlobNumbering(pTag, int EntType);

};
