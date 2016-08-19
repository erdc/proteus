#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <gmi_null.h>
#include <ma.h>
#include <maShape.h>
#include <apfMDS.h>
#include <PCU.h>
#include <SimUtil.h>
#include <SimModel.h>

#include <iostream>
#include <fstream>

#include "MeshAdaptPUMI.h"

MeshAdaptPUMIDrvr::MeshAdaptPUMIDrvr(double Hmax, double Hmin, int NumIter,
    const char* sfConfig, const char* maType,const char* logType, double targetError, double targetElementCount)
{
  m = 0;
  PCU_Comm_Init();
  PCU_Protect();
  Sim_readLicenseFile(0);
  SimModel_start();
  hmin=Hmin; hmax=Hmax;
  numIter=NumIter;
  nAdapt=0;
  nEstimate=0;
  if(PCU_Comm_Self()==0)
     printf("MeshAdapt: Setting hmax=%lf, hmin=%lf, numIters(meshadapt)=%d\n",
       hmax, hmin, numIter);
  global[0] = global[1] = global[2] = global[3] = 0;
  local[0] = local[1] = local[2] = local[3] = 0;
  size_iso = 0;
  size_scale = 0;
  size_frame = 0;
  err_reg = 0;
  errRho_reg = 0;
  errRel_reg = 0;
  gmi_register_mesh();
  gmi_register_sim();
  gmi_register_null();
  approximation_order = 2;
  integration_order = 3;//approximation_order * 2;
  total_error = 0.0;
  errRho_max = 0.0;
  rel_err_total = 0.0;
  exteriorGlobaltoLocalElementBoundariesArray = NULL;
  size_field_config = sfConfig;
  modelFileName = NULL; 
  adapt_type_config = maType;
  logging_config = logType;
  has_gBC = false;
  target_error = targetError;
  target_element_count = targetElementCount;
  domainVolume = 0.0;
}

MeshAdaptPUMIDrvr::~MeshAdaptPUMIDrvr()
{
  freeField(err_reg);
  freeField(errRho_reg);
  freeField(errRel_reg);
  freeField(size_iso);
  freeField(size_scale);
  freeField(size_frame);
  SimModel_stop();
  Sim_unregisterAllKeys();
}


static bool ends_with(std::string const& str, std::string const& ext)
{
  return str.size() >= ext.size() &&
         str.compare(str.size() - ext.size(), ext.size(), ext) == 0;
}

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

#include <MeshSim.h>
#include <SimMeshTools.h>
#define FACE 2
pAManager SModel_attManager(pModel model);
/*
Temporary function used to read in BC from Simmetrix Model
*/
int MeshAdaptPUMIDrvr::getSimmetrixBC()
{
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
  return 0;
} 

int MeshAdaptPUMIDrvr::willAdapt() //THRESHOLD needs to be defined
{
  double THRESHOLD = 0;//target_error;
  int adaptFlag=0;
  int assertFlag;

  if(total_error > THRESHOLD){
    adaptFlag = 1;
  }
  assertFlag = adaptFlag;
  PCU_Add_Ints(&assertFlag,1);
  assert(assertFlag ==0 || assertFlag == PCU_Proc_Peers());
  return adaptFlag;
}

int MeshAdaptPUMIDrvr::adaptPUMIMesh()
{
  if (size_field_config == "farhad")
      calculateAnisoSizeField();
  else if (size_field_config == "alvin"){
      removeBCData();
      double t1 = PCU_Time();
      getERMSizeField(total_error);
      freeField(err_reg); //mAdapt will throw error if not destroyed. what about free?
      freeField(errRho_reg); //mAdapt will throw error if not destroyed. what about free?
      freeField(errRel_reg); //mAdapt will throw error if not destroyed. what about free?
      double t2 = PCU_Time();
    if(comm_rank==0 && logging_config == "on"){
      std::ofstream myfile;
      myfile.open("error_estimator_timing.txt", std::ios::app );
      myfile << t2-t1<<std::endl;
      myfile.close();
    }
  
  }  
  else if (size_field_config == "isotropic")
    testIsotropicSizeField();
  else {
    std::cerr << "unknown size field config " << size_field_config << '\n';
    abort();
  }
  //m->destroyTag(fluxtag[1]); m->destroyTag(fluxtag[2]); m->destroyTag(fluxtag[3]);
  delete [] exteriorGlobaltoLocalElementBoundariesArray;
  exteriorGlobaltoLocalElementBoundariesArray = NULL;

  for (int d = 0; d <= m->getDimension(); ++d)
    freeNumbering(local[d]);
  /// Adapt the mesh
  ma::Input* in;
  if(adapt_type_config=="anisotropic")
    in = ma::configure(m, size_scale, size_frame);
  else
    in = ma::configure(m, size_iso);
  ma::validateInput(in);
  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->maximumImbalance = 1.05;
  in->maximumIterations = numIter;
  in->shouldSnap = false;
  in->shouldFixShape = true;
  double mass_before = getTotalMass();

  double t1 = PCU_Time();
  ma::adapt(in);
  double t2 = PCU_Time();

  freeField(size_iso);
  freeField(size_frame);
  freeField(size_scale);
  m->verify();
  double mass_after = getTotalMass();
  PCU_Add_Doubles(&mass_before,1);
  PCU_Add_Doubles(&mass_after,1);
  if(comm_rank==0 && logging_config=="on"){
/*
    std::ios::fmtflags saved(std::cout.flags());
    std::cout<<std::setprecision(15)<<"Mass Before "<<mass_before<<" After "<<mass_after<<" diff "<<mass_after-mass_before<<std::endl;
    std::cout.flags(saved);
*/
    std::ofstream myfile;
    myfile.open("adapt_timing.txt", std::ios::app);
    myfile << t2-t1<<std::endl;
    myfile.close();

    std::ofstream mymass;
    mymass.open("mass_check.txt", std::ios::app);
    mymass <<std::setprecision(15)<<mass_before<<","<<mass_after<<","<<mass_after-mass_before<<std::endl;
    mymass.close();
  }
  if(size_field_config=="alvin"){
    if (has_gBC)
      getSimmetrixBC();
    if(logging_config=="on"){
      char namebuffer[20];
      sprintf(namebuffer,"pumi_postadapt_%i",nAdapt);
      apf::writeVtkFiles(namebuffer, m);
    }
  }
  nAdapt++; //counter for number of adapt steps
  return 0;
}

double MeshAdaptPUMIDrvr::getMinimumQuality()
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
  return mass;
}
