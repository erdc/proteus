#include "mesh.h"
#include "cmeshToolsModule.h"
#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>

class MeshAdaptPUMIDrvr{
 
  public:
  MeshAdaptPUMIDrvr(double, double, int, const char*, const char*,const char*); 
  ~MeshAdaptPUMIDrvr();

  int loadModelAndMesh(const char* modelFile, const char* meshFile);

  //Functions to construct proteus mesh data structures
  int constructFromSerialPUMIMesh(Mesh& mesh);
  int constructFromParallelPUMIMesh(Mesh& mesh, Mesh& subdomain_mesh);
 
  int updateMaterialArrays(Mesh& mesh, int bdryID, int GeomTag);

  int transferFieldToPUMI(const char* name, double const* inArray, int nVar, int nN);
  int transferFieldToProteus(const char* name, double* outArray, int nVar, int nN);
  int transferPropertiesToPUMI(double* rho_p, double* nu_p,double* g_p);
  int transferBCtagsToProteus(int* tagArray, int idx, int* ebN, int* eN_global, double* fluxBC);
  int transferBCsToProteus();
  int commuSizeField();
  int adaptPUMIMesh();

  void numberLocally();
  int localNumber(apf::MeshEntity* e);
  int dumpMesh(Mesh& mesh);

  int calculateSizeField();
  int calculateAnisoSizeField();
  int testIsotropicSizeField();
  int getERMSizeField(double err_total,double rel_err_total);
  double getMinimumQuality();
  double getTotalMass();
  double getMPvalue(double field_val,double val_0, double val_1);
  apf::Field* getViscosityField(apf::Field* voff);

  double hmax, hmin;
  int numIter;
  int nAdapt; //counter for number of adapt steps
  int nsd; //number of spatial dimensions

  //Element Residual Method
  void get_local_error();
  void computeDiffusiveFlux(apf::Mesh*m,apf::Field* voff, apf::Field* visc,apf::Field* pref, apf::Field* velf);
  //void getBoundaryFlux(apf::Mesh* m, apf::MeshEntity* ent, apf::Field* voff, apf::Field* visc,apf::Field* pref, apf::Field* velf, double * endflux);
  void getBoundaryFlux(apf::Mesh* m, apf::MeshEntity* ent, double * endflux);
  int getSimmetrixBC();
  void simmetrixBCreloaded(const char* modelFile);
  void removeBCData();
  char* modelFileName; 
  
  //tags used to identify types of BC
  apf::MeshTag* BCtag;
  apf::MeshTag* DBCtag[4];
  apf::MeshTag* fluxtag[4];

  int *exteriorGlobaltoLocalElementBoundariesArray;

  //Approximation/Integration order
  int approximation_order; //what order polynomial (hierarchic is 2nd order)
  int integration_order; //determines number of integration points
  int num_quadrature; 
  int num_quarature_boundary;

  private: 
  apf::Mesh2* m;
  int comm_size, comm_rank;

  double rho[2], nu[2];
  double g[3];
  apf::MeshTag* diffFlux;
  apf::GlobalNumbering* global[4];
  apf::Numbering* local[4];
  apf::Field* err_reg; //error field from ERM
  apf::Field* rel_err; //relative error field
  /* this field stores isotropic size */
  apf::Field* size_iso;
  /* these fields store anisotropic size */
  apf::Field* size_scale;
  apf::Field* size_frame;

  int constructGlobalNumbering(Mesh& mesh);
  int constructGlobalStructures(Mesh& mesh);

  int constructElements(Mesh& mesh);
  int constructNodes(Mesh& mesh);
  int constructBoundaries(Mesh& mesh);
  int constructEdges(Mesh& mesh);
  int constructMaterialArrays(Mesh& mesh);

  void freeField(apf::Field*& f);
  void freeNumbering(apf::Numbering*& n);

  static void averageToEntity(apf::Field* ef, apf::Field* vf,
      apf::MeshEntity* ent);

  std::string size_field_config;
  std::string adapt_type_config;
  std::string logging_config;
};


