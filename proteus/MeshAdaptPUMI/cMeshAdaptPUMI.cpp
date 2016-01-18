#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <ma.h>
#include <maShape.h>
#include <apfMDS.h>
#include <PCU.h>
#include <SimUtil.h>
#include <SimModel.h>

#include "MeshAdaptPUMI.h"

MeshAdaptPUMIDrvr::MeshAdaptPUMIDrvr(double Hmax, double Hmin, int NumIter,
    const char* sfConfig, const char* maType)
{
  m = 0;
  PCU_Comm_Init();
  PCU_Protect();
  Sim_readLicenseFile(0);
  SimModel_start();
  hmin=Hmin; hmax=Hmax;
  numIter=NumIter;
  nAdapt=0;
  if(PCU_Comm_Self()==0)
     printf("MeshAdapt: Setting hmax=%lf, hmin=%lf, numIters(meshadapt)=%d\n",
       hmax, hmin, numIter);
  global[0] = global[1] = global[2] = global[3] = 0;
  local[0] = local[1] = local[2] = local[3] = 0;
  size_iso = 0;
  size_scale = 0;
  size_frame = 0;
  err_reg = 0;
  gmi_register_mesh();
  gmi_register_sim();
  approximation_order = 2;
  integration_order = 3;//approximation_order * 2;
  casenum = 2;
  exteriorGlobaltoLocalElementBoundariesArray = NULL;
  size_field_config = sfConfig;
  geomFileName = NULL; 
  modelFileName = NULL; 
  meshFileName = NULL; 
  adapt_type_config = maType;
}

MeshAdaptPUMIDrvr::~MeshAdaptPUMIDrvr()
{
  freeField(size_iso);
  freeField(size_scale);
  freeField(size_frame);
  SimModel_stop();
  Sim_unregisterAllKeys();
}

int MeshAdaptPUMIDrvr::loadModelAndMesh(const char* modelFile, const char* meshFile)
{
  m = apf::loadMdsMesh(modelFile, meshFile);
  m->verify();
  comm_size = PCU_Comm_Peers();
  comm_rank = PCU_Comm_Self();
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
     std::cout<<"Found case, setting the model"<<std::endl;
     AttCase_setModel(acase,model);
  } else {
      std::cout<<"Case not found, rename case to geom\n"<<std::endl;
      exit(1);
  }
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
  int bcFlag[nsd+1] = {0,1,1,1};

  //assign a label to the BC type tag
  char label[4][9],labelflux[6],type_flag;
  for(int idx=0;idx<4;idx++)
  {
    if(idx == 0) sprintf(&type_flag,"p");
    else if(idx == 1) sprintf(&type_flag,"u");
    else if(idx == 2) sprintf(&type_flag,"v");
    else if(idx == 3) sprintf(&type_flag,"w");
    sprintf(label[idx],"BCtype_%c",type_flag);
    BCtag[idx] = m->createIntTag(label[idx],1);
std::cout<<"Boundary label "<<label[idx]<<std::endl;
    if(idx>0) sprintf(labelflux,"%c_flux",type_flag);
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
        fluxtag[idx]= m->createDoubleTag(labelflux,numqpt);
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
          for(int idx=1;idx<nsd+1;idx++)
          {
            m->setIntTag(fEnt,BCtag[idx],&bcFlag[idx]);
            m->setDoubleTag(fEnt,fluxtag[idx],data[idx]); //set the quadrature points
          }
          apf::destroyMeshElement(testElem);
        } //end if on model
        else
        {
          for(int idx=1;idx<nsd+1;idx++)
          {
            int dummy = 0;
            m->setIntTag(fEnt,BCtag[idx],&dummy);
          }
        }
      }//end loop over attributes
      if(nF==0)
      {
          for(int idx=1;idx<nsd+1;idx++)
          {
            int dummy = 0;
            m->setIntTag(fEnt,BCtag[idx],&dummy);
          }
      }
    } 
  }//end while
  m->end(fIter);
  AMAN_release( attmngr );
  std::cout<<"Finished reading and storing diffusive flux BCs\n"; 
  return 0;
} 

void MeshAdaptPUMIDrvr::simmetrixBCreloaded(const char* modelFile)
{ 
  if(modelFileName == NULL)
  {
    modelFileName=(char *) malloc(sizeof(char) * strlen(modelFile));
    strcpy(modelFileName,modelFile);
  }
  for(int i=0;i<comm_size;i++){
    if(comm_rank==i)
      getSimmetrixBC();
    PCU_Barrier();
  }
}

int MeshAdaptPUMIDrvr::adaptPUMIMesh()
{
  if (size_field_config == "farhad")
    calculateAnisoSizeField();
  else if (size_field_config == "alvin")
    get_local_error();
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
  ma::Input* in = ma::configure(m, size_scale, size_frame);
  ma::validateInput(in);
  in->shouldRunPreParma = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->maximumIterations = numIter;
  in->shouldSnap = false;
  in->shouldFixShape = true;
  double mass_before = getTotalMass();
  ma::adapt(in);
  freeField(size_frame);
  freeField(size_scale);
  m->verify();
  double mass_after = getTotalMass();
  std::ios::fmtflags saved(std::cout.flags());
  std::cout<<std::setprecision(15)<<"Before "<<mass_before<<" After "<<mass_after<<" diff "<<mass_after-mass_before<<std::endl;
  std::cout.flags(saved);
/*
  if(size_field_config=="alvin")
    simmetrixBCreloaded(modelFileName);
*/
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
