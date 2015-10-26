#include "mesh.h"
#include "cmeshToolsModule.h"
#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>

class MeshAdaptPUMIDrvr{
 
  public:
  MeshAdaptPUMIDrvr(double, double, int); 
  ~MeshAdaptPUMIDrvr();

  int loadModelAndMesh(const char* modelFile, const char* meshFile);

  //Functions to construct proteus mesh data structures
  int ConstructFromSerialPUMIMesh(Mesh& mesh);
  int ConstructFromParallelPUMIMesh(Mesh& mesh, Mesh& subdomain_mesh);
 
  int UpdateMaterialArrays(Mesh& mesh, int bdryID, int GeomTag);

  //Fields
  //Transfer Boundary Conditions
  int TransferSolutionToPUMI(double* inArray, int nVar, int nN);
  int TransferSolutionToProteus(double* outArray, int nVar, int nN);
  int TransferPropertiesToPUMI(double* rho_p, double* nu_p);
  int TransferBCtagsToProteus();
  int TransferBCsToProteus();
  int CommuSizeField();
  int AdaptPUMIMesh();

  void numberLocally();
  int localNumber(apf::MeshEntity* e);
  int dumpMesh(Mesh& mesh);

  int CalculateSizeField();
  int CalculateAnisoSizeField();
  int getERMSizeField(double err_total);

  double hmax, hmin;
  int numIter;
  int nAdapt; //counter for number of adapt steps

  //Element Residual Method
  void get_local_error();
  //tags used to identify types of BC
  apf::MeshTag* BCtag[4];
  apf::MeshTag* DBCtag[4];
  apf::MeshTag* fluxtag[4];

  //Approximation/Integration order
  int approximation_order; //what order polynomial (hierarchic is 2nd order)
  int integration_order; //determines number of integration points

  private: 
  apf::Mesh2* m;
  int comm_size, comm_rank;
  int elms_owned, faces_owned, edges_owned, vtx_owned;
  int numVar;

  double rho[2], nu[2];
  apf::GlobalNumbering* global[4];
  apf::Numbering* local[4];
  apf::Field* solution;
  apf::Field* fluxBC;
  apf::Field* DBC;
  apf::Field* err_reg; //error field from ERM
  /* this field stores isotropic size */
  apf::Field* size_iso;
  /* these fields store anisotropic size */
  apf::Field* size_scale;
  apf::Field* size_frame;

  int ConstructGlobalNumbering(Mesh& mesh);
  int ConstructGlobalStructures(Mesh& mesh);

  int ConstructElements(Mesh& mesh);
  int ConstructNodes(Mesh& mesh);
  int ConstructBoundaries(Mesh& mesh);
  int ConstructEdges(Mesh& mesh);
  int ConstructMaterialArrays(Mesh& mesh);

  void freeField(apf::Field*& f);
  void freeNumbering(apf::Numbering*& n);
};


