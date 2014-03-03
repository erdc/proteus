#include <cstdlib>
//#include "MA.h"
#include "pumi_mesh.h"
#include "pumi.h"
#include "pumi_geom.h"
#include "mesh.h"
#include "cmeshToolsModule.h"
#include "MeshTools.h"

void TransferTopSCOREC(pPList oldEnts, pPList newRgn, void *userData, modType mtype, pEntity ent);

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

  //Fields
  int TransferSolutionToPUMI(double* inArray, int nVar, int nN);
  int TransferSolutionToProteus(double* outArray, int nVar, int nN);
  int CalculateSizeField();
  int CommuSizeField();
  int AdaptPUMIMesh();
  int MeshAdaptPUMI();

  int InterpolateSolutionE( pEdge edge, double xi[2], int field_size, pTag pTagTag, double* result);
  int TransferBottomE(pPList parent, pPList fresh, pPList VtxstoHandle, modType mtype);

  private: 
  pMeshMdl PUMI_MeshInstance;
  pGeomMdl PUMI_GModel;
  pPart PUMI_Part;
  std::vector<pPart> PUMI_Parts;
  int comm_size, comm_rank;
  int elms_owned, faces_owned, edges_owned, vtx_owned;
  int numVar;

  pTag GlobNumberTag;
  pTag SolutionTag, SFTag;
  std::vector<std::vector<double> > SolutionVec;

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
  int DeleteMeshEntIDs();
};
