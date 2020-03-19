#include <gmi.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <ma.h>
#include <maShape.h>
#include <apfMDS.h>
#include <PCU.h>
#include <apf.h>

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>

#include "MeshAdaptPUMI.h"
#include <sam.h>
#include <samSz.h>

#include <apfShape.h>

#include "PyEmbeddedFunctions.h"

extern double dt_err;
#ifdef PROTEUS_USE_SIMMETRIX
//PROTEUS_USE_SIMMETRIX is a compiler macro that indicates whether Simmetrix libraries are used
//This is defined in proteus/config/default.py and is contingent on the existence of a SIM_INCLUDE_DIR path
  #include <gmi_sim.h>
  #include <SimUtil.h>
  #include <SimModel.h>
  #include <MeshSim.h>
  #include <SimMeshTools.h>
  #define FACE 2
  pAManager SModel_attManager(pModel model);
#endif

/** 
 * \file cMeshAdaptPUMI.cpp
 * \ingroup MeshAdaptPUMI 
 @{ 
*/
//MeshAdaptPUMIDrvr::MeshAdaptPUMIDrvr(double Hmax, double Hmin, double HPhi,int AdaptMesh, int NumIter, int NumAdaptSteps,const char* sfConfig, const char* maType,const char* logType, double targetError, double targetElementCount,int reconstructedFlag,double maxAspectRatio, double gradingFact)
MeshAdaptPUMIDrvr::MeshAdaptPUMIDrvr()
/**
 * MeshAdaptPUMIDrvr is the highest level class that handles the interface between Proteus and the PUMI libraries
 * See MeshAdaptPUMI.h for the list of class variables/functions/objects
 * This is the constructor for the class
*/
{
  m = 0;
  PCU_Comm_Init();
  PCU_Protect();

#ifdef PROTEUS_USE_SIMMETRIX
  Sim_readLicenseFile(0);
  SimModel_start();
  gmi_register_sim();
#endif
  nAdapt=0;
  nTriggers=0;
  nEstimate=0;
  global[0] = global[1] = global[2] = global[3] = 0;
  local[0] = local[1] = local[2] = local[3] = 0;
  size_iso = 0;
  size_scale = 0;
  size_frame = 0;
  err_reg = 0;
  vmsErrH1 = 0;
  errRho_reg = 0;
  errRel_reg = 0;
  error_reference = 0;
  gmi_register_mesh();
  gmi_register_null();
  approximation_order = 2;
  integration_order = 3;
  total_error = 0.0;
  errRho_max = 0.0;
  rel_err_total = 0.0;
  exteriorGlobaltoLocalElementBoundariesArray = NULL;
  modelFileName = NULL; 
  adapt_type_config = "test"; //refers to isotropic or anisotropic, feature effectively turned off for now
  has_gBC = false;
  target_element_count = 0;//targetElementCount;
  domainVolume = 0.0;
  THRESHOLD = 0.0;
  initialReconstructed = 0;
  maxAspect = 2.0;//maxAspectRatio;
  hasIBM=hasInterface=hasVMS=hasERM=hasAniso=hasAnalyticSphere=useProteus=useProteusAniso=0;
}

MeshAdaptPUMIDrvr::~MeshAdaptPUMIDrvr()
/**
 * Destructor for MeshAdaptPUMIDrvr
 */
{
/*
  freeField(err_reg);
  freeField(vmsErrH1);
  freeField(errRho_reg);
  freeField(errRel_reg);
  freeField(error_reference);
  freeField(size_iso);
  freeField(size_scale);
  freeField(size_frame);
*/
/*
  if(isReconstructed){
    free(modelVertexMaterial);
    free(modelBoundaryMaterial);
    free(modelRegionMaterial);
    //m->destroyNative();
    //gmi_destroy(m->getModel());
    //apf::destroyMesh(m);
  }
*/
  PCU_Comm_Free();
#ifdef PROTEUS_USE_SIMMETRIX
  SimModel_stop();
  Sim_unregisterAllKeys();
#endif
}


static bool ends_with(std::string const& str, std::string const& ext)
{
  return str.size() >= ext.size() &&
         str.compare(str.size() - ext.size(), ext.size(), ext) == 0;
}
/**
 * @brief Load the mesh and model for SCOREC libraries
 *
 * The default filetypes are .dmg (model) and .smb (mesh), but can support GMSH meshes (.msh) and Simmetrix models (.smd) and meshes (.sms)
 * GMSH models are not used and .null filenames and passed instead.
 * Each of the the GMSH and Simmetrix filetypes can be converted into a SCOREC filetype via tools in scorec
 * Diffusive flux boundary conditions are supported with Simmetrix models and can be passed into the error estimator 
 */

int MeshAdaptPUMIDrvr::loadModelAndMesh(const char* modelFile, const char* meshFile)
{
  comm_size = PCU_Comm_Peers();
  comm_rank = PCU_Comm_Self();
  if (ends_with(meshFile, ".msh")){
    m = apf::loadMdsFromGmsh(gmi_load(modelFile), meshFile);
    std::cout<<"Boundary Condition functionality has not been built in for gmsh yet.\n";
  }
  else if (ends_with(modelFile,".smd")){
    m = apf::loadMdsMesh(modelFile, meshFile);
    modelFileName=(char *) malloc(sizeof(char) * strlen(modelFile));
    strcpy(modelFileName,modelFile);
    getSimmetrixBC();
  }
  else{
    m = apf::loadMdsMesh(modelFile, meshFile);
  }

  m->verify();
  return 0;
}

int MeshAdaptPUMIDrvr::loadMeshForAnalytic(const char* meshFile,double* boxDim,double* sphereCenter, double radius)
{
  //assume analytic 
  comm_size = PCU_Comm_Peers();
  comm_rank = PCU_Comm_Self();
  m = apf::loadMdsMesh(".null", meshFile);
  m->verify();

  //create analytic geometry 
  gmi_model* testModel = createSphereInBox(boxDim,sphereCenter,radius);
  m->verify();


/*
  apf::writeVtkFiles("afterAnalytic",m);
  std::cout<<"test Model "<<testModel<<" mesh model "<<m->getModel()<<std::endl;
  std::abort();
*/
  return 0;
}



int MeshAdaptPUMIDrvr::getSimmetrixBC()
/**
 * @brief Function used to read in diffusive flux BC from Simmetrix Model.
 *
 * Simmetrix BCs set via the GUI are read-in through the Simmetrix API.
 * The values are stored as apf tags on the mesh.
 * BCs are not currently supported with any other type of model.
 */
{

#ifdef PROTEUS_USE_SIMMETRIX
  pGModel model = 0;
  model=GM_load(modelFileName,NULL,NULL);

  pAManager attmngr = SModel_attManager(model);
  pACase acase = AMAN_findCaseByType(attmngr, "problem definition");
  if (acase){
     if(comm_rank==0)std::cout<<"Found case, setting the model"<<std::endl;
     AttCase_setModel(acase,model);
     has_gBC=true;
  AttCase_associate(acase,NULL);

  pGFace gFace;
  GFIter gfIter = GM_faceIter(model);
  pAttribute Att[GM_numFaces(model)];
  int attMap[GM_numFaces(model)];
  int nF=0;
  
  char strAtt[2][25] = {"traction vector","comp3"};
  int modelEntTag;

  while(gFace = GFIter_next(gfIter))
  {
    if(GEN_attrib((pGEntity)gFace,strAtt[0]))
    { 
      modelEntTag=GEN_tag((pGEntity)gFace);
      Att[nF]=GEN_attrib((pGEntity)gFace,strAtt[0]);
      attMap[nF] = modelEntTag;
      nF++;
    }
  }
  GFIter_delete(gfIter);
  
  apf::MeshIterator* fIter = m->begin(FACE);
  apf::MeshEntity* fEnt;
  apf::Vector3 evalPt;
  int numqpt=0;
  const int nsd = 3;
  const int bcFlag[nsd+1] = {0,1,1,1};

  //assign a label to the BC type tag
  char label[9],labelflux[4][9],type_flag;
  
  sprintf(label,"BCtype");
  BCtag = m->createIntTag(label,4);
  for(int idx=0;idx<4;idx++)
  {
    if(idx == 0) sprintf(&type_flag,"p");
    else if(idx == 1) sprintf(&type_flag,"u");
    else if(idx == 2) sprintf(&type_flag,"v");
    else if(idx == 3) sprintf(&type_flag,"w");
    if(idx>0) sprintf(labelflux[idx],"%c_flux",type_flag);
  }

  while(fEnt = m->iterate(fIter))
  {
    apf::ModelEntity* me=m->toModel(fEnt);
    modelEntTag = m->getModelTag(me);
    apf::ModelEntity* boundary_face = m->findModelEntity(FACE,modelEntTag);
    if(numqpt==0)
    {
      apf::MeshElement* testElem = apf::createMeshElement(m,fEnt);
      numqpt = apf::countIntPoints(testElem,integration_order);
      for(int idx=1;idx<nsd+1;idx++)
        fluxtag[idx]= m->createDoubleTag(labelflux[idx],numqpt);
      apf::destroyMeshElement(testElem);
    }
    if(me==boundary_face)
    {
      for(int i=0;i<nF;i++)
      {
        if(attMap[i]==modelEntTag)
        {
          apf::MeshElement* testElem = apf::createMeshElement(m,fEnt);
          double data[nsd+1][numqpt];
          for(int k=0; k<numqpt;k++)
          {
            apf::getIntPoint(testElem,integration_order,k,evalPt);
            apf::Vector3 evalPtGlobal;
            apf::mapLocalToGlobal(testElem,evalPt,evalPtGlobal);
            double evalPtSim[nsd];
            evalPtGlobal.toArray(evalPtSim);
            for(int j=0;j<nsd;j++)
              data[j+1][k]=AttributeTensor1_evalDS((pAttributeTensor1)Att[i], j,evalPtSim);
          }
          m->setIntTag(fEnt,BCtag,&(bcFlag[0]));
          for(int idx=1;idx<nsd+1;idx++)
          {
            m->setDoubleTag(fEnt,fluxtag[idx],data[idx]); //set the quadrature points
          }
          apf::destroyMeshElement(testElem);
          break;
        } //end if on model
        else
        {
          int dummy[4] = {0,0,0,0};
          m->setIntTag(fEnt,BCtag,&(dummy[0]));
        }
      }//end loop over attributes
      if(nF==0)
      {
          int dummy[4] = {0,0,0,0};
          m->setIntTag(fEnt,BCtag,&(dummy[0]));
      }
    } 
  }//end while
  m->end(fIter);
  AMAN_release( attmngr );
  } else {
      if(comm_rank==0)
        std::cout<<"Case not found, no BCs?\n"<<std::endl;
      //exit(1);
  }

  if(comm_rank==0)std::cout<<"Finished reading and storing diffusive flux BCs\n"; 
#endif
  return 0;
} 
static int countTotal(apf::Mesh* m, int dim)
{
  int total = apf::countOwned(m, dim);
  PCU_Add_Ints(&total, 1);
  return total;
}


int MeshAdaptPUMIDrvr::willErrorAdapt() 
/**
 * @brief Looks at the estimated error and determines if mesh adaptation is necessary.
 *
 * Function used to define whether a mesh needs to be adapted based on the error estimator
 * The return value is a flag indicating whether the mesh will not (0) or will be (1) adapted 
 * The THRESHOLD will be set to the error estimate after the wind-up step, but is currently 0
 * Assertion is set to ensure that all ranks in a parallel execution will enter the adapt stage
 */
{
  int adaptFlag=0;
  int assertFlag;

  //get current size field
  apf::Field* currentField;

/*
  if(!size_iso) //if no previous size field
  {
    currentField = samSz::isoSize(m);
  }
  else //otherwise use previous size field, remember to reflect this in interfaceAdapt or collapse to a single function
  {
    currentField  = apf::createFieldOn(m, "currentField", apf::SCALAR);
    apf::copyData(currentField,size_iso);
  }
*/
  //currentField  = apf::createFieldOn(m, "currentField", apf::SCALAR);
  //apf::copyData(currentField,size_iso);
  currentField = samSz::isoSize(m);

  //get error-based size field
  getERMSizeField(total_error);
  apf::Field* errorField = sizeFieldList.front();
  sizeFieldList.pop(); //remove this size field from the queue


  apf::Field *errorTriggered = apf::createLagrangeField(m, "errorTriggered", apf::SCALAR, 1);
  //determine if desired mesh is contained in current mesh
  apf::MeshEntity* ent;
  apf::MeshIterator* it = m->begin(0);
  while( (ent = m->iterate(it)) )
  {
    double h_current = apf::getScalar(currentField,ent,0);
    double h_needed = apf::getScalar(errorField,ent,0);
    if(h_current>h_needed*1.5){
      adaptFlag+=1;        
      apf::setScalar(errorTriggered,ent,0,h_current/h_needed*1.0);
      //apf::writeVtkFiles("willErrorAdapt", m);
      //std::cout<<"What is the ent? "<<localNumber(ent)<<std::endl;
      //std::exit(1);
      //break;
    }
    else
      apf::setScalar(errorTriggered,ent,0,-1);

  }//end while
  m->end(it);

  //modify error field to be the ratio
  it = m->begin(0);
  while( (ent = m->iterate(it)) )
  {
    double h_current = apf::getScalar(currentField,ent,0);
    double h_needed = apf::getScalar(errorField,ent,0);
    apf::setScalar(errorField,ent,0,h_current/h_needed*1.0);
  }//end while
  m->end(it);

  assertFlag = adaptFlag;
  PCU_Add_Ints(&assertFlag,1);
  //assert(assertFlag ==0 || assertFlag == PCU_Proc_Peers());

  if(assertFlag>0)
  {
    double totalNodes = countTotal(m,0);
    double triggeredPercentage = assertFlag*100.0/totalNodes;
    char buffer[50];
    sprintf(buffer,"Need to error adapt %f%%",triggeredPercentage);
    logEvent(buffer,3);

/*
    if(nTriggers%10==0)
    {
        char namebuffer[50];
        sprintf(namebuffer,"needErrorAdapt_%i",nTriggers);
        apf::writeVtkFiles(namebuffer, m);
    }
*/
    nTriggers++;
  }

  apf::destroyField(currentField);
  //apf::destroyField(errorField);
  apf::destroyField(errorTriggered);

  return assertFlag;
}

int MeshAdaptPUMIDrvr::willErrorAdapt_reference() 
{
  int adaptFlag=0;
  int assertFlag;

  getERMSizeField(total_error);
  sizeFieldList.pop(); //remove this size field from the queue

  //if(m->findField("errorTriggered"))
  //  apf::destroyField(m->findField("errorTriggered"));
  apf::Field* errorTriggered;
  if(m->findField("errorTriggered"))
    errorTriggered = m->findField("errorTriggered");
  else
    errorTriggered = apf::createField(m,"errorTriggered",apf::SCALAR,apf::getVoronoiShape(m->getDimension(),1));
 
  apf::Field* error_current = m->findField("VMSH1");
  apf::Field* error_reference=NULL;


  if(m->findField("errorRate"))
    apf::destroyField(m->findField("errorRate"));
  apf::Field* errorRateField = apf::createField(m,"errorRate",apf::SCALAR,apf::getVoronoiShape(m->getDimension(),1));

  //need to set the error trigger field to 0
  apf::MeshEntity* ent;
  apf::MeshIterator* it;

  if(nTriggers == 0)
  {
    it = m->begin(m->getDimension());
    while( (ent = m->iterate(it) ) )
    {
        apf::setScalar(errorTriggered,ent,0,-1.0);
        apf::setScalar(errorRateField,ent,0,0.0);
    }
    m->end(it);
    logEvent("SET ERROR TRIGGERED!",4);
  }


  //need to have size field here just to define it

  //set the error reference field to be the current error field for next time otherwise, retrieve the error_reference field
  if(!m->findField("error_reference"))
  {  
    T_reference = T_current;
    error_reference = apf::createField(m,"error_reference",apf::SCALAR,apf::getVoronoiShape(nsd,1));
    apf::copyData(error_reference,error_current);
    logEvent("SUCCESSFULLY COPIED!",4);

    return 0;
  }
  else
  {
    error_reference = m->findField("error_reference");
  }

  double dt_step = dt_err; //global variable imported from VMS

  if(m->findField("sizeRatio"))
    apf::destroyField(m->findField("sizeRatio"));
  apf::Field* sizeRatioField = apf::createField(m,"sizeRatio",apf::SCALAR,apf::getVoronoiShape(m->getDimension(),1));

  it = m->begin(m->getDimension());

  while( (ent = m->iterate(it) ) )
  {   
    double err_local_current = apf::getScalar(error_current,ent,0);
    double err_local_ref = apf::getScalar(error_reference,ent,0);
    double errorRate = (err_local_current-err_local_ref)/(T_current-T_reference);
    apf::setScalar(errorRateField,ent,0,errorRate);
    double err_predict = errorRate*delta_T_next + err_local_current;

    //double sizeRatio_local = apf::getScalar(m->findField("sizeRatio"),ent,0);
    //if(errorRate > 0 && (err_predict > target_error) && sizeRatio_local>2.0 && nTriggers>=5)
    //double sizeRatio_local_predict = pow(err_predict/target_error,2.0/(2*1+2));

    double h_old;
    double h_new;
    apf::MeshElement* element = apf::createMeshElement(m, ent);

    if (m->getDimension() == 2)
      h_old = apf::computeLargestHeightInTri(m,ent);
    else
      h_old = apf::computeShortestHeightInTet(m,ent);
    h_new = h_old * pow((target_error / err_predict),2.0/(2.0*(1.0)+nsd));
    //clamp h_new and then compare against h_old
    if(h_new > hmax)
        h_new = hmax;
    if(h_new < hmin)
        h_new = hmin;
    double sizeRatio_local_predict = h_old/h_new;
    apf::setScalar(sizeRatioField,ent,0,sizeRatio_local_predict);
    apf::destroyMeshElement(element);

    //apf::setScalar(sizeRatioField,ent,0,size_ratio_local_predict);
    if(errorRate > 0 && sizeRatio_local_predict>2.0 && nTriggers>= (0.1*numAdaptSteps))
    {
        if(apf::getScalar(errorTriggered,ent,0)==1)
        {
            adaptFlag = 1;
            apf::setScalar(errorTriggered,ent,0,2.0);   
        }
        else
            apf::setScalar(errorTriggered,ent,0,1.0);
    }
    else
        apf::setScalar(errorTriggered,ent,0,-1.0);
  }
  m->end(it);

  //upper bound on adapt steps
  if(nTriggers>=2*numAdaptSteps)
  {
    logEvent("Adapting because upper limit was reached.",4);
    adaptFlag = 1;
  }

  assertFlag = adaptFlag;
  PCU_Add_Ints(&assertFlag,1);
  nTriggers++;

  //if adapt, modify the error field to be predictive
  //otherwise, store the current error field to be the next reference field
  if(assertFlag > 0)
  {
    logEvent("Need to error based adapt!!!",4);
/*
    it = m->begin(m->getDimension());
    while( (ent = m->iterate(it) ) )
    {   
        double err_local_current = apf::getScalar(error_current,ent,0);
        double err_local_ref = apf::getScalar(error_reference,ent,0);
        //double errorRate = (err_local_current-err_local_ref)/(T_current-T_reference);
        double errorRate = apf::getScalar(errorRateField,ent,0);
    
        //double err_predict = errorRate*dt_step*numAdaptSteps + err_local_current;
        double err_predict = errorRate*dt_step*numAdaptSteps/10.0 + err_local_current;
        if(err_predict > err_local_current)
        {
            apf::setScalar(error_current,ent,0,err_predict);
        }
    }
    m->end(it);
*/
  }
  else
  {
    it = m->begin(m->getDimension());
    while( (ent = m->iterate(it) ) )
    {   
        double err_local_current = apf::getScalar(error_current,ent,0);
        //set the reference field for the next step
        apf::setScalar(error_reference,ent,0,err_local_current);
    }

  }
  
  T_reference = T_current;

  return assertFlag;
}



int MeshAdaptPUMIDrvr::willAdapt() 
//Master function that calls other adapt-trigger functions
{
  int adaptFlag = 0;
  //TODO: need to separate IBM from interface to allow for two-phase with IBM
  if(hasInterface or hasIBM)
  {
    adaptFlag += willInterfaceAdapt(); 
    std::cout<<"willInterfaceAdapt "<<adaptFlag<<std::endl;
  }
  if(hasVMS)
  {
    adaptFlag += willErrorAdapt_reference();
    std::cout<<"willErrorAdapt "<<adaptFlag<<std::endl;
  }
  //TODO: need to do AllReduce here instead of in each function
  if(adaptFlag > 0)
    adaptFlag = 1;

  if(adaptFlag == 0)
  {
    //allocated in transfer fields... this is not a good way of doing things, but don't know how to pass a numpy array without having to allocate memory just yet
    free(rho);
    free(nu);
  }

  return adaptFlag;
}


int MeshAdaptPUMIDrvr::willInterfaceAdapt() 
//Does banded adapt need to happen for an isotropic mesh?
//I need to loop over all mesh edges and determine if the edge intersects the blending region.
//If so, need to check the size values on the edge-adjacent vertices.
//If either size value is greater than h_interface*1.5, then we know we need to adapt
{
  int adaptFlag=0;
  int assertFlag;

  //get current size field
  apf::Field* currentField;
/*
  if(!size_iso) //if no previous size field
  {
    currentField = samSz::isoSize(m);
  }
  else //otherwise use previous size field, remember to reflect this in interfaceAdapt or collapse to a single function
  {
    currentField  = apf::createFieldOn(m, "currentField", apf::SCALAR);
    apf::copyData(currentField,size_iso);
  }
*/

  currentField = samSz::isoSize(m);
  double edgeRatio = 1.5; //need to be taken from MeshAdapt library

  //get banded size field
  double L_band = (N_interface_band)*hPhi;
  calculateSizeField(L_band);
  apf::Field* interfaceField = sizeFieldList.front();
  sizeFieldList.pop(); //destroy this size field


  //determine if desired mesh is contained in current mesh
  apf::MeshEntity* ent;
  apf::MeshIterator* it = m->begin(0);
  while( (ent = m->iterate(it)) )
  {
    double h_current = apf::getScalar(currentField,ent,0);
    double h_needed = apf::getScalar(interfaceField,ent,0);
    if(h_current/h_needed > edgeRatio){
      adaptFlag=1;        
      break;
    }  
  }//end while

  assertFlag = adaptFlag;
  PCU_Add_Ints(&assertFlag,1);
  //assert(assertFlag ==0 || assertFlag == PCU_Proc_Peers());

  apf::destroyField(currentField);
  apf::destroyField(interfaceField);

  return assertFlag;
}



int MeshAdaptPUMIDrvr::adaptPUMIMesh(const char* inputString)
/**
 * @brief Function used to trigger adaptation
 *
 * Inputs are the type of size-field that is desired:
 * "interface" refers to explicitly adapting to the current interface position based on the level-set field
 * "ERM" refers to the error-residual method and using error estimates to determine how a mesh should be adapted
 *  Within ERM, there is an isotropic and anisotropic configuration. The anisotropic configuration requires more development.
 *  "Isotropic" refers to a uniform refinement based on a global hmin size and is primarily used for testing
 *  Predictive load balancers Zoltan and ParMA are used in combination for before, during, and after adapt with a preset tolerance of imbalance
 *  The nAdapt counter is iterated to track how many times a mesh is adapted
 */
{
  if (hasAniso)
      calculateAnisoSizeField();
  if(hasERM)
  {
      assert(err_reg);
      removeBCData();
      double t1 = PCU_Time();
      getERMSizeField(total_error);
      double t2 = PCU_Time();
      if(comm_rank==0 && logging_config == "on"){
        std::ofstream myfile;
        myfile.open("error_estimator_timing.txt", std::ios::app );
        myfile << t2-t1<<std::endl;
        myfile.close();
      }
  } 
  if((hasVMS) && std::string(inputString)==""){
    assert(vmsErrH1);
    getERMSizeField(total_error);
  }
  if(hasAnalyticSphere)
  {
    setSphereSizeField();
  }
  if(hasInterface || hasIBM)
  {
    double L_band = (N_interface_band+1)*hPhi;
    calculateSizeField(L_band);
    if(nAdapt>1 && std::string(inputString)=="interface")
        predictiveInterfacePropagation();
    else if(nAdapt > 2 && std::string(inputString)=="")
        predictiveInterfacePropagation();
  }
  if (useProteus)
    size_iso = m->findField("proteus_size");
  if (useProteusAniso){
    size_frame = m->findField("proteus_sizeFrame");
    size_scale = m->findField("proteus_sizeScale");
    adapt_type_config = "anisotropic";
  }
/*
  if (hasTest)
    testIsotropicSizeField();
*/
/*
  if(size_field_config == "uniform"){
      //special situation where I only care about err_reg
      freeField(errRho_reg); 
      freeField(errRel_reg); 
  }
*/
  std::cout<<"before problem\n";
  isotropicIntersect();
  std::cout<<"after problem\n";

  if(logging_config=="on"){
    char namebuffer[50];
    sprintf(namebuffer,"pumi_preadapt_%i",nAdapt);
    apf::writeVtkFiles(namebuffer, m);
    //sprintf(namebuffer,"beforeAnisotropicAdapt%i_.smb",nAdapt);
    //m->writeNative(namebuffer);
  }

  if(hasERM){
      //MeshAdapt error will be thrown if region fields are not freed
      freeField(err_reg); 
      freeField(errRho_reg); 
      freeField(errRel_reg); 
  }
  if(hasVMS){
    freeField(vmsErrH1);
    if(PCU_Comm_Self()==0) std::cout<<"cleared VMS field\n";
  }

  // These are relics from an attempt to pass BCs from proteus into the error estimator.
  // They maybe useful in the future.
  //m->destroyTag(fluxtag[1]); m->destroyTag(fluxtag[2]); m->destroyTag(fluxtag[3]);
  delete [] exteriorGlobaltoLocalElementBoundariesArray;
  exteriorGlobaltoLocalElementBoundariesArray = NULL;

  for (int d = 0; d <= m->getDimension(); ++d)
    freeNumbering(local[d]);

  apf::Field* adaptSize;
  apf::Field* adaptFrame;

  /// Adapt the mesh
  ma::Input* in;
  if(size_field_config == "uniform"){
    in = ma::configureUniformRefine(m);
    in->shouldFixShape=false;
  }
  else{
    assert(size_iso || (size_scale && size_frame));
    if(adapt_type_config=="anisotropic" || size_field_config== "interface"){
     //in = ma::configure(m, size_scale, size_frame);
      adaptSize  = apf::createFieldOn(m, "adapt_size", apf::VECTOR);
      adaptFrame = apf::createFieldOn(m, "adapt_frame", apf::MATRIX);
      apf::copyData(adaptSize, size_scale);
      apf::copyData(adaptFrame, size_frame);
      in = ma::configure(m, adaptSize, adaptFrame);
    }
    else{
      adaptSize  = apf::createFieldOn(m, "adapt_size", apf::SCALAR);
      apf::copyData(adaptSize, size_iso);
      in = ma::configure(m, adaptSize);
    }
  }

  ma::validateInput(in);
  in->shouldRunPreZoltan = true;
  in->shouldRunMidZoltan = true;
  in->shouldRunPostZoltan = true;
  in->maximumImbalance = 1.05;
  in->maximumIterations = numIter;
  if(size_field_config == "meshQuality")
  {
    in->shouldSnap = true;
    in->shouldTransferParametric=true;
  }
  else
    in->shouldSnap = false;
  //in->goodQuality = 0.16;//0.027;
  //double mass_before = getTotalMass();
  
  double t1 = PCU_Time();
  //ma::adapt(in);
  ma::adaptVerbose(in);
  double t2 = PCU_Time();

  m->verify();
  //double mass_after = getTotalMass();
  //PCU_Add_Doubles(&mass_before,1);
  //PCU_Add_Doubles(&mass_after,1);
  if(comm_rank==0 && logging_config=="on"){
/*
    std::ios::fmtflags saved(std::cout.flags());
    std::cout<<std::setprecision(15)<<"Mass Before "<<mass_before<<" After "<<mass_after<<" diff "<<mass_after-mass_before<<std::endl;
    std::cout.flags(saved);
    std::ofstream myfile;
    myfile.open("adapt_timing.txt", std::ios::app);
    myfile << t2-t1<<std::endl;
    myfile.close();
    std::ofstream mymass;
    mymass.open("mass_check.txt", std::ios::app);
    mymass <<std::setprecision(15)<<mass_before<<","<<mass_after<<","<<mass_after-mass_before<<std::endl;
    mymass.close();
*/
  }

  if(hasERM){
    if (has_gBC)
      getSimmetrixBC();
  }
  if(logging_config=="on"){
    char namebuffer[50];
    sprintf(namebuffer,"pumi_postadapt_%i",nAdapt);
    apf::writeVtkFiles(namebuffer, m);
    //sprintf(namebuffer,"afterAnisotropicAdapt%i_.smb",nAdapt);
    //m->writeNative(namebuffer);
  }
  //isReconstructed = 0; //this is needed to maintain consistency with the post-adapt conversion back to Proteus
  apf::destroyField(adaptSize);
  if(adapt_type_config=="anisotropic")
    apf::destroyField(adaptFrame);
  nAdapt++; //counter for number of adapt steps

  if(logging_config=="debugRestart")
    m->writeNative("DEBUG_restart.smb");

  nTriggers=0;
  return 0;
}

double MeshAdaptPUMIDrvr::getMinimumQuality()
/**
 * @brief Function used to get the worst element quality in the mesh.
 *
 * Measures the quality via SCOREC library; returns the minimum quality
 * Meant to be used to trigger adaptation, but has not been implemented yet
 */
{
  ma::SizeField* isf = new ma::IdentitySizeField(m);
  apf::MeshIterator* it = m->begin(m->getDimension());
  apf::MeshEntity* e;
  double minq = 1;
  while ((e = m->iterate(it)))
    minq = std::min(minq, ma::measureElementQuality(m, isf, e));
  m->end(it);
  delete isf;
  return PCU_Min_Double(minq);
}

double MeshAdaptPUMIDrvr::getTotalMass()
/**
 * @brief Function to track total mass of the domain.
 *
 * Returns total mass across all ranks
 */
{
  apf::Field* voff = m->findField("vof");
  assert(voff);
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(m->getDimension());
  double mass = 0.0;
  while ((e = m->iterate(it))) {
    apf::MeshElement* elem = apf::createMeshElement(m,e);
    apf::Element* voff_elem = apf::createElement(voff, elem);
    int int_order = 4;
    for(int l = 0; l < apf::countIntPoints(elem, int_order); ++l) {
      apf::Vector3 qpt;
      apf::getIntPoint(elem,int_order,l,qpt);
      double vof_val = apf::getScalar(voff_elem,qpt);
      double rho_val = getMPvalue(vof_val,rho[0],rho[1]);
      double weight = apf::getIntWeight(elem,int_order,l);
      apf::Matrix3x3 J;
      apf::getJacobian(elem,qpt,J); //evaluate the Jacobian at the quadrature point
      double Jdet = apf::getJacobianDeterminant(J,m->getDimension());
      mass += rho_val*weight*Jdet;
    }
    apf::destroyElement(voff_elem);
    apf::destroyMeshElement(elem);
  }
  m->end(it);
  PCU_Add_Doubles(&mass,1);
  return mass;
}
/** @} */

//Save mesh with solution

void MeshAdaptPUMIDrvr::writeMesh(const char* meshFile)
{
  m->writeNative(meshFile);
  //apf::writeVtkFiles(meshFile,m);
}

//Clean mesh of all fields and tags

void MeshAdaptPUMIDrvr::cleanMesh()
{
  //destroy all fields...
    
  for(int i =0;i<m->countFields();i++)
  {
    apf::Field* sample = m->getField(i);
    freeField(sample);
  }
  //std::cout<<"find field " <<m->getField(m->countFields())<<" how many fields? "<<m->countFields()<<std::endl;
  //std::cout<<"is velocity_old here? "<<m->findField("velocity_old")<<std::endl;
  apf::Field* sample = m->findField("velocity_old");
  freeField(sample);

  sample = m->findField("vof_old");
  freeField(sample);
  sample = m->findField("ls_old");
  freeField(sample);
  sample = m->findField("phi");
  freeField(sample);
  sample = m->findField("phi_old");
  freeField(sample);
  sample = m->findField("phi_old_old");
  freeField(sample);
  sample = m->findField("phid_old");
  freeField(sample);
  sample = m->findField("phiCorr");
  freeField(sample);
  sample = m->findField("phiCorr_old");
  freeField(sample);
  sample = m->findField("phiCorr_old_old");
  freeField(sample);
  sample = m->findField("p_old");
  freeField(sample);
  sample = m->findField("p");
  freeField(sample);
  sample = m->findField("p_old_old");
  freeField(sample);

  sample = m->findField("VMSH1");
  freeField(sample);
  sample = m->findField("VMSL2");
  freeField(sample);

  //destroy all tags
  apf::DynamicArray<apf::MeshTag*> listTags;
  m->getTags(listTags);
  int nTags = listTags.getSize();
  int numDim = m->getDimension();
  for(int i=0; i < nTags; i++)
  {
    std::string ignoreString ("proteus_number");
    std::string tagName (m->getTagName(listTags[i]));
    if(tagName.find(ignoreString) != std::string::npos)
    {
      //do nothing
    }
    else{
      for(int j=0;j<(numDim+1);j++)
      {
        apf::MeshIterator* it = m->begin(j);
        apf::MeshEntity* ent;
        while( (ent = m->iterate(it)) )
        {
          if(m->hasTag(ent,listTags[i]))
            m->removeTag(ent,listTags[i]);
        }
        m->end(it);
      }
      m->destroyTag(listTags[i]);
    }
  }
}

void MeshAdaptPUMIDrvr::set_nAdapt(int numberAdapt)
{
  nAdapt = numberAdapt;
  return;
}

int MeshAdaptPUMIDrvr::setAdaptProperties(std::vector<std::string> sizeInputs,bool in_adapt, double in_hmax,double in_hmin,double in_hphi, int in_numAdaptSteps, double in_targetError, double in_gradingFactor, bool in_logging, int in_numIterations)
{
    for (std::vector<std::string>::iterator it = sizeInputs.begin() ; it != sizeInputs.end(); ++it)
    {
        if(*it == "ibm")
            hasIBM=1;
        if(*it == "interface")
            hasInterface=1;
        if(*it == "error_vms")
            hasVMS=1;
        if(*it == "error_erm")
            hasERM=1;
    }
    hmin=in_hmin; 
    hmax=in_hmax; 
    hPhi=in_hphi;
    numIter=in_numIterations;
    adaptMesh = in_adapt;
    numAdaptSteps = in_numAdaptSteps;
    if(PCU_Comm_Self()==0)
        printf("MeshAdapt: Setting hmax=%lf, hmin=%lf, numIters(meshadapt)=%d\n",
         hmax, hmin, numIter);
    logging_config = in_logging;
    target_error = in_targetError;
    gradingFactor = in_gradingFactor;
   
    return 0;
}
