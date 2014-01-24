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
  int deleteProteusMesh(Mesh& mesh);

  int readGeomModel(const std::string &acis_geom_file_name);
  int readPUMIMesh(const char* SMS_fileName);
  int helloworld(const char* hello) { std::cout << hello << "\n";}

  //Functions to construct proteus mesh data structures
//  int ConstructFromPUMIMesh(Mesh& mesh);
  int ConstructFromPUMIMesh(CMesh *cmesh);
  int ConstructElements(Mesh& mesh);
  int ConstructNodes(Mesh& mesh);
  int ConstructBoundaries(Mesh& mesh);
  int ConstructEdges(Mesh& mesh);
  int ConstructMaterialArrays(Mesh& mesh);

  int MeshAdaptPUMI();

  private: 
  pMeshMdl PUMI_MeshInstance;
  pGeomMdl PUMI_GModel;
  std::vector<pPart> PUMI_Parts;
  pPart PUMI_Part;
//  pMMeshAdaptPUMI MA_Drvr;

};
