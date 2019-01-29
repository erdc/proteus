#include "mesh.h"
#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>

/**
   \file MeshAdaptPUMI.h
   \defgroup MeshAdaptPUMI MeshAdaptPUMI
   \brief A C-library of functions that uses SCOREC tools for error estimation and mesh adaptation
   @{
*/

/**
  @brief Class that handles the interface between SCOREC tools and Proteus
*/

class MeshAdaptPUMIDrvr{
 
  public:
  MeshAdaptPUMIDrvr(double, double, double, int, int, int, const char*, const char*,const char*,double,double,int,double,double); 
  ~MeshAdaptPUMIDrvr();

  int loadModelAndMesh(const char* modelFile, const char* meshFile); //load the model and mesh
  void writeMesh(const char* meshFile);

  //Functions to construct proteus mesh data structures
  int reconstructFromProteus(Mesh& mesh, Mesh& globalMesh,int hasModel);
  int reconstructFromProteus2(Mesh& mesh, int* isModelVert, int* bFaces);
  int constructFromSerialPUMIMesh(Mesh& mesh);
  int constructFromParallelPUMIMesh(Mesh& mesh, Mesh& subdomain_mesh);
  int updateMaterialArrays(Mesh& mesh,int dim, int bdryID, int GeomTag);
  int updateMaterialArrays(Mesh& mesh);
  int updateMaterialArrays2(Mesh& mesh);
  void numberLocally();
  int localNumber(apf::MeshEntity* e);
  int dumpMesh(Mesh& mesh);

  //Functions used to transfer proteus model data structures
  int transferModelInfo(int*numGeomEntities,int*edges,int*faces,int*mVertex2Model,int*mEdgeVertex2Model,int*mBoundary2Model,int nMaxSegments);
  int numSegments;
  int* edgeList;
  int* faceList;
  int* meshVertex2Model, *meshEdge2Model, *meshBoundary2Model;
  int numModelEntities[4];

  //Functions used to transfer information between PUMI and proteus
  int transferFieldToPUMI(const char* name, double const* inArray, int nVar, int nN);
  int transferFieldToProteus(const char* name, double* outArray, int nVar, int nN);
  int transferPropertiesToPUMI(double* rho_p, double* nu_p,double* g_p, double deltaT, double interfaceBandSize);
  //int transferBCtagsToProteus(int* tagArray, int idx, int* ebN, int* eN_global, double* fluxBC);
  //int transferBCsToProteus();

  //MeshAdapt functions
  int willAdapt();
  int adaptPUMIMesh();
  int calculateSizeField();
  int calculateAnisoSizeField();
  int testIsotropicSizeField();
  int getERMSizeField(double err_total);
  int gradeMesh();

  //Quality Check Functions
  double getMinimumQuality();
  double getTotalMass();

  //Functions that help facilitate computations
  double getMPvalue(double field_val,double val_0, double val_1); //get the multiphase value of physical properties
  apf::Field* getViscosityField(apf::Field* voff); //derive a field of viscosity based on VOF field


  //Public Variables
  double hmax, hmin, hPhi; //bounds on mesh size
  int numIter; //number of iterations for MeshAdapt
  int nAdapt; //counter for number of adapt steps
  int nEstimate; //counter for number of error estimator calls
  int nsd; //number of spatial dimensions
  int maxAspect; //maximum aspect ratio
  int adaptMesh; //is adaptivity turned on?
  int numAdaptSteps; //Number adaptivity
  double N_interface_band; //number of elements in half-band around interface
  double gradingFactor;

  //User Inputs
  std::string size_field_config; //What type of size field: interface, ERM, isotropic
  std::string adapt_type_config; //What type of adapt for ERM: isotropic or anisotropic
  std::string logging_config; // Logging on or off

  //Element Residual Method
  void get_local_error(double& total_error);
  void computeDiffusiveFlux(apf::Mesh*m,apf::Field* voff, apf::Field* visc,apf::Field* pref, apf::Field* velf);
  void getBoundaryFlux(apf::Mesh* m, apf::MeshEntity* ent, double * endflux);
  int getSimmetrixBC();
  void removeBCData();
  char* modelFileName; 
  
  //VMS
  void get_VMS_error(double& total_error);
  
  //tags used to identify types of BC
  apf::MeshTag* BCtag;
  apf::MeshTag* DBCtag[4];
  apf::MeshTag* fluxtag[4];

  int *exteriorGlobaltoLocalElementBoundariesArray;

  //Approximation/Integration order
  int approximation_order; //what order polynomial (hierarchic is 2nd order)
  int integration_order; //determines number of integration points
  int num_quadrature; 
  int num_quadrature_boundary;
  double total_error;
  double errRho_max;
  double rel_err_total;

  //Mesh Reconstruction
  int isReconstructed;
  int initialReconstructed;
  int* modelVertexMaterial; 
  int* modelBoundaryMaterial;
  int* modelRegionMaterial;
  int numModelOffsets[4];
  int numModelTotals[4];

  private: 
  apf::Mesh2* m;
  int comm_size, comm_rank;

  double rho[2], nu[2];
  double g[3];
  double delta_T;
  apf::MeshTag* diffFlux;
  apf::GlobalNumbering* global[4];
  apf::Numbering* local[4];
  apf::Field* err_reg; //error field from ERM
  apf::Field* vmsErrH1; //error field for VMS
  apf::Field* errRho_reg; //error-density field from ERM
  apf::Field* errRel_reg; //relative error field from ERM
  /* this field stores isotropic size */
  apf::Field* size_iso;
  /* these fields store anisotropic size and metric tensor */
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
  void volumeAverageToEntity(apf::Field* ef, apf::Field* vf,
      apf::MeshEntity* ent);

  bool has_gBC; //boolean for having global boundary conditions
  double target_error; //computed from get_local_error()
  int target_element_count; //How many elements in the mesh are desired?
  double domainVolume; //Volume of the domain
  double THRESHOLD; //threshold of error before adapt
};


/** @} */
