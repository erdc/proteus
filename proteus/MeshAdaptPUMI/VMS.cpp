#include "MeshAdaptPUMI.h"
#include <PCU.h>
#include <petscksp.h>

#include <apf.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <apfDynamicMatrix.h>
#include <apfNumbering.h>

#include <iostream>
#include <fstream>

//Global variables used to make it easier to pass these variables from MeshAdaptPUMIDrvr
extern int approx_order; //shape function order
extern int int_order; //integration order
extern double nu_0,nu_1,rho_0,rho_1;
double dt_err;

struct Inputs{
  apf::Vector3 vel_vect;  
  apf::Matrix3x3 gij;
  apf::Matrix3x3 grad_vel;
  apf::Vector3 grad_pres;
  double visc_val;
  double density;
  int nsd;
  apf::MeshElement* element;
  apf::Element* pres_elem;
  apf::Element* visc_elem;
  apf::Element* velo_elem;
  apf::Element* vof_elem;
  apf::Element* velo_elem_old;
  apf::Matrix3x3 KJ;
  double* g;
};

double get_nu_err(struct Inputs info);
apf::Vector3 getResidual(apf::Vector3 qpt,struct Inputs &info);
apf::Vector3 getResidual_debug(apf::Vector3 qpt,struct Inputs &info);


inline void getProps(double*rho,double*nu,double deltaT)
//Function used to transfer MeshAdaptPUMIDrvr variables into global variables
{
  rho_0 = rho[0];
  nu_0 = nu[0];      
  rho_1 = rho[1];
  nu_1 = nu[1];
  dt_err = deltaT;
  std::cout<<"THis is delta T right now! "<<deltaT<<" and rho and nu "<<rho[0]<<" "<<nu[0]<<std::endl;
  return;
}

apf::Matrix3x3 getKJ(int nsd)
{
//Jacobian correction matrix for proper metric tensor
    if(nsd==2)
    {
        apf::Matrix3x3 KJtemp(1.0,0.5,0.0,0.5,1.0,0.0,0.0,0.0,1.0); //Jacobian correction matrix for proper metric tensor
        return KJtemp;
    }
    else
    {
        double c1 = pow(2,1./3.); //1.259921049894873e+00;
        double c2 = c1/2.0; //6.299605249474365e-01;
        apf::Matrix3x3 KJtemp(c1,c2,c2,c2,c1,c2,c2,c2,c1);
        return KJtemp;
    }
}

double MeshAdaptPUMIDrvr::getTotalLinearMomentum()
{
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(m->getDimension());
  apf::Field* velf = m->findField("velocity");
  double x_momentum = 0.0;
  double y_momentum = 0.0;
  while ((e = m->iterate(it))) {
    apf::MeshElement* elem = apf::createMeshElement(m,e);
    apf::Element* vel_elem = apf::createElement(velf, elem);
    int int_order = 4;
    double density = rho[localNumber(e)];
    for(int l = 0; l < apf::countIntPoints(elem, int_order); ++l) {
      apf::Vector3 qpt;
      apf::getIntPoint(elem,int_order,l,qpt);
      apf::Vector3 velocity;
      apf::getVector(vel_elem,qpt,velocity);
      double weight = apf::getIntWeight(elem,int_order,l);
      apf::Matrix3x3 J;
      apf::getJacobian(elem,qpt,J); //evaluate the Jacobian at the quadrature point
      double Jdet = apf::getJacobianDeterminant(J,m->getDimension());
      //momentum += density*weight*Jdet;
      x_momentum += velocity[0]*density*weight*Jdet;
      y_momentum += velocity[1]*density*weight*Jdet;
    }
    apf::destroyElement(vel_elem);
    apf::destroyMeshElement(elem);
  }
  m->end(it);
  //PCU_Add_Doubles(&mass,1);
  std::cout<<nAdapt<<" What is the momentum? "<<x_momentum<<" "<<y_momentum<<std::endl;
  return x_momentum;
}



void MeshAdaptPUMIDrvr::get_VMS_error(double &total_error_out) 
{
  if(PCU_Comm_Self()==0)
    std::cout<<"The beginning of the VMS\n";
  getProps(rho,nu,delta_T);
  approx_order = approximation_order; 
  int_order = integration_order;
  nsd = m->getDimension();

  //***** Get Solution Fields First *****//
  apf::Field* voff = m->findField("vof");
  assert(voff);
  apf::Field* velf = m->findField("velocity");
  assert(velf);
  apf::Field* pref = m->findField("p");
  assert(pref);
  apf::Field* velf_old;
  if(m->findField("velocity_old")!=NULL)
    velf_old = m->findField("velocity_old");
  else{
    velf_old = velf;
    if(PCU_Comm_Self()==0)
      std::cout<<"WARNING: old velocity field not found. Will proceed as if unsteady term is 0. \n";
    dt_err = 1.0;
  }
  assert(velf_old);
  //*****               *****//
  if(PCU_Comm_Self()==0)
    std::cout<<"Got the solution fields\n";
  //***** Compute the viscosity field *****//
  apf::Field* visc = getViscosityField(voff);
  if(PCU_Comm_Self()==0)
    std::cout<<"Got viscosity fields \n";
  freeField(vmsErrH1);

  apf::Field* vmsErr = apf::createField(m,"VMSL2",apf::SCALAR,apf::getVoronoiShape(nsd,1));
  //vmsErrH1 = apf::createField(m,"VMSH1",apf::SCALAR,apf::getVoronoiShape(nsd,1));
  vmsErrH1 = apf::createField(m,"VMSH1",apf::SCALAR,apf::getVoronoiShape(nsd,1));
  if(m->findField("nu_err"))
    apf::destroyField(m->findField("nu_err"));
  nuErr  = apf::createField(m,"nu_err",apf::SCALAR,apf::getVoronoiShape(nsd,1));
  if(m->findField("strongResidual"))
    apf::destroyField(m->findField("strongResidual"));
  strongResidual  = apf::createField(m,"strongResidual",apf::VECTOR,apf::getVoronoiShape(nsd,1));
  
  if(PCU_Comm_Self()==0)
    std::cout<<"Created the error fields\n";
  //Start computing element quantities
  int numqpt; //number of quadrature points
  int nshl; //number of local shape functions
  int elem_type; //what type of topology
  double weight; //value container for the weight at each qpt
  double Jdet;
  //apf::EntityShape* elem_shape;
  apf::Vector3 qpt; //container for quadrature points
  apf::MeshElement* element;
  apf::Element* visc_elem, *pres_elem,*velo_elem,*vof_elem,*velo_elem_old;
  apf::Matrix3x3 J; //actual Jacobian matrix
  apf::Matrix3x3 invJ; //inverse of Jacobian
  apf::Matrix3x3 KJ = getKJ(nsd);
  apf::Matrix3x3 gij; //actual Jacobian matrix


  //apf::DynamicMatrix invJ_copy;
  //apf::NewArray <apf::DynamicVector> shdrv;
  //apf::NewArray <apf::DynamicVector> shgval_copy;
  
  apf::MeshIterator* iter = m->begin(nsd); //loop over elements
  apf::MeshEntity* ent;

  int count = 0 ; //for debugging purposes
  //Loop over elements and compute the VMS error in the L_2 norm
  double VMSerrTotalL2 = 0.0;
  double VMSerrTotalH1 = 0.0;
  while( (ent = m->iterate(iter)) ){

    element = apf::createMeshElement(m,ent);
    pres_elem = apf::createElement(pref,element);
    velo_elem = apf::createElement(velf,element);
    velo_elem_old = apf::createElement(velf_old,element);
    visc_elem = apf::createElement(visc,element);
    vof_elem = apf::createElement(voff,element);
    
    double strongResidualTauL2;
    double tau_m;
    double areaCheck=0.0;
    numqpt=apf::countIntPoints(element,int_order);

    //Start the quadrature loop
    for(int k=0;k<numqpt;k++){
      apf::getIntPoint(element,int_order,k,qpt); //get a quadrature point and store in qpt
      apf::getJacobian(element,qpt,J); //evaluate the Jacobian at the quadrature point
      weight = apf::getIntWeight(element,int_order,k);
      J = apf::transpose(J); //Is PUMI still defined in this way?
      if(nsd==2)
        J[2][2] = 1.0; //this is necessary to avoid singular matrix
      invJ = invert(J);
      Jdet=fabs(apf::getJacobianDeterminant(J,nsd)); 
      gij = apf::transpose(invJ)*(KJ*invJ);
      
/*
    if(count == 10 && k==0){
      apf::Adjacent adj;
      m->getAdjacent(ent,0,adj);
      for(int i=0;i<adj.getSize();i++){
        apf::Vector3 pt;
        m->getPoint(adj[i],0,pt);
        std::cout<<std::setprecision(15)<<std::scientific;
        std::cout<<"The point is "<<pt[0]<<" "<<pt[1]<<" "<<pt[2]<<std::endl;
      }
      apf::Matrix3x3 testJ;
      apf::Vector3 pt1,pt2,pt3;
      m->getPoint(adj[0],0,pt1);
      m->getPoint(adj[1],0,pt2);
      m->getPoint(adj[2],0,pt3);
      testJ[0][0] = pt2[0]-pt1[0]; 
      testJ[0][1] = pt3[0]-pt1[0]; 
      testJ[1][0] = pt2[1]-pt1[1]; 
      testJ[1][1] = pt3[1]-pt1[1]; 
      testJ[2][2] = 1;
     
      std::cout<<"Jacobian "<<J<<" inverse "<<invJ<<" testJ "<<testJ<<std::endl;
      std::cout<<"Jdet "<<Jdet<<" area "<<apf::measure(m,ent)<<std::endl;
      for(int i=0;i<adj.getSize();i++){
        double presCheck = apf::getScalar(pref,adj[i],0);
        std::cout<<"Pressure is "<<presCheck<<std::endl;
      }
      apf::Vector3 gradCheck;
      apf::getGrad(pres_elem,qpt,gradCheck);
      std::cout<<"grad pres is "<<gradCheck<<std::endl;

      for(int i=0;i<adj.getSize();i++){
        apf::Vector3 velCheck; 
        apf::getVector(velf,adj[i],0,velCheck);
        std::cout<<"velocity is "<<velCheck<<std::endl;
      }
      apf::Matrix3x3 gradVelCheck;
      apf::getVectorGrad(velo_elem,qpt,gradVelCheck);
      std::cout<<"grad vel is "<<gradVelCheck<<std::endl;
    }
*/
      
      double pressure = apf::getScalar(pres_elem,qpt);
      apf::Vector3 grad_pres;
      apf::getGrad(pres_elem,qpt,grad_pres);
      double visc_val = apf::getScalar(visc_elem,qpt);
      apf::Vector3 vel_vect;
      apf::getVector(velo_elem,qpt,vel_vect);
      apf::Matrix3x3 grad_vel;
      apf::getVectorGrad(velo_elem,qpt,grad_vel);
      grad_vel = apf::transpose(grad_vel);

      double C1 = 2.0; //constants for stabilization term
      double C2 = 36.0; 

      double stabTerm1 = 0.0;
      double stabTerm2 = 0.0;
      for(int i=0;i<nsd;i++){
        for(int j=0;j<nsd;j++){
          stabTerm1 += vel_vect[i]*gij[i][j]*vel_vect[j];
          stabTerm2 += gij[i][j]*gij[i][j];
        }
      }
      stabTerm2 = C2*visc_val*visc_val*stabTerm2;
      tau_m = 1/sqrt(stabTerm1 + stabTerm2);

/*      
      if(count == 0){
        //std::cout<<"This is the tau "<<tau_m<<" velocity "<<vel_vect<<" grad vel "<<grad_vel<<" metric "<<gij<<" numdim "<<nsd<<" viscosity "<<visc_val<<" "<<KJ<<" gradP"<< grad_pres<<" stabTerms "<<stabTerm1<<" "<<stabTerm2<<std::endl;
        areaCheck += weight*Jdet;
      }
*/

      //compute residual
      apf::Vector3 tempConv;
      apf::Vector3 tempDiff;
      tempConv.zero();
      tempDiff.zero();
      for(int i=0;i<nsd;i++){
        for(int j=0;j<nsd;j++){
          tempConv[i] = tempConv[i] + vel_vect[j]*grad_vel[i][j];
        }
      }
      double density = getMPvalue(apf::getScalar(vof_elem,qpt),rho_0,rho_1);
      apf::Vector3 tempResidual = (tempConv + grad_pres/density);
      double tempVal = tempResidual.getLength();
      strongResidualTauL2 = tau_m*tau_m*tempVal*tempVal*weight*Jdet;
/*
      if(count == 10 && k==0){
        std::cout<<" tempConv "<<tempConv<<" density "<<density<<" tempResidual "<<tempResidual<<" tempVal "<<tempVal<<std::endl;
      }
*/

    } //end qpt loop

    double VMSerrL2 = sqrt(strongResidualTauL2);
    apf::setScalar(vmsErr,ent,0,VMSerrL2);
/*
    if(count == 10){
      std::cout<<"This is the error "<< sqrt(strongResidualTauL2)<<std::endl;
      std::cout<<"This is the area check "<<areaCheck<<" true area is "<<apf::measure(m,ent);
    }
*/
    
    //H1 error compute nu_err at centroid and compute residual
    qpt[0] = 0.0; qpt[1] = 0.0; qpt[2] = 0.0; //ensure it is initialized to 0
    for(int i=0; i<nsd; i++)
      qpt[i] = 1./(nsd+1.0);

    //double density = getMPvalue(apf::getScalar(vof_elem,qpt),rho_0,rho_1);
    //double density = getMPvalue(apf::getScalar(vof_elem,qpt),998.2,1.205);
    //std::cout<<" double check density "<<rho[localNumber(ent)]<<" "<<density<<std::endl;
    double density = rho[localNumber(ent)];
    //std::abort();
    Inputs info;
    info.element = element;
    info.pres_elem = pres_elem;
    info.visc_elem = visc_elem;
    info.velo_elem = velo_elem;
    info.vof_elem = vof_elem;
    info.velo_elem_old = velo_elem_old;
    info.KJ = KJ;
    info.nsd = nsd;
    info.density = density;
    info.g = &g[0];

/*
    double pressure = apf::getScalar(pres_elem,qpt);
    apf::Vector3 grad_pres;
    apf::getGrad(pres_elem,qpt,grad_pres);
    double visc_val = apf::getScalar(visc_elem,qpt);

    apf::Vector3 vel_vect;
    apf::getVector(velo_elem,qpt,vel_vect);
    apf::Matrix3x3 grad_vel;
    apf::getVectorGrad(velo_elem,qpt,grad_vel);
    grad_vel = apf::transpose(grad_vel);
      //compute residual
      apf::Vector3 tempConv;
      apf::Vector3 tempDiff;
      tempConv.zero();
      tempDiff.zero();
      for(int i=0;i<nsd;i++){
        for(int j=0;j<nsd;j++){
          tempConv[i] = tempConv[i] + vel_vect[j]*grad_vel[i][j];
        }
      }
      double density = getMPvalue(apf::getScalar(vof_elem,qpt),rho_0,rho_1);
      apf::Vector3 tempResidual = (tempConv + grad_pres/density);
      double tempVal = tempResidual.getLength();
*/
    info.visc_val = nu[localNumber(ent)];
    apf::Vector3 tempResidual = getResidual(qpt,info);
    double tempVal = tempResidual.getLength();
    if(nAdapt == 4 && (localNumber(ent) == 800 || localNumber(ent)==490))
    {
        std::cout<<nAdapt<<" element number "<<localNumber(ent)<<std::endl;
        getResidual_debug(qpt,info);
    }
    if(nAdapt == 5 && (localNumber(ent) == 6404 || localNumber(ent) == 6554))
    {
        std::cout<<nAdapt<<" element number "<<localNumber(ent)<<std::endl;
        getResidual_debug(qpt,info);
    }
    apf::getJacobian(element,qpt,J);
    if(nsd==2)
      J[2][2] = 1.0; //this is necessary to avoid singular matrix
    invJ = invert(J);
    Jdet=fabs(apf::getJacobianDeterminant(J,nsd)); 
    gij = apf::transpose(invJ)*(KJ*invJ);
    info.gij = gij;

    double nu_err = get_nu_err(info);
    double VMSerrH1 = nu_err*tempVal*sqrt(apf::measure(m,ent));
    apf::setScalar(nuErr,ent,0,nu_err);
    apf::setVector(strongResidual,ent,0,tempResidual);

/*
    if(nAdapt==5)
    {
        if(localNumber(ent) == 1061 || localNumber(ent) == 2475)
        {
            std::cout<<"Element number "<<localNumber(ent)<<" has strong residual "<<tempResidual<<" and nu_err "<<nu_err<<std::endl;
            std::abort();
        }
    }
*/

    apf::Vector3 pt_centroid = apf::getLinearCentroid(m,ent);
    
    if(fabs(pt_centroid[0]-0.508)<0.01 && fabs(pt_centroid[1]-0.01556)<0.01 && PCU_Comm_Self()==2 && nAdapt>6)
    {
        std::cout<<"error vms is "<<VMSerrH1<<std::endl;
        std::cout<<"What is the residual "<<tempResidual<<" Jacobian "<<J<<" gij "<<gij<< " density "<<density<<" nu "<< info.visc_val<<std::endl;
        std::cout<<"vel "<<info.vel_vect<<" grad_vel "<<info.grad_vel<<" grad pres "<<info.grad_pres<<std::endl;
        std::cout<<"What is the area of triangle? "<<apf::measure(m,ent)<<" dt_err "<<dt_err<<std::endl;
        std::cout<<"nu err "<<nu_err<<std::endl;
        apf::Vector3 vel_vect_old;
        apf::getVector(info.velo_elem_old,qpt,vel_vect_old);

        std::cout<<"vel_vect_old "<< vel_vect_old<<std::endl;
        apf::Adjacent adjVerts;
        m->getAdjacent(ent,0,adjVerts);
        apf::Vector3 pt_coord;
        for(int i =0; i<3; i++)
        {
            m->getPoint(adjVerts[i],0,pt_coord);
            std::cout<<"coordinates "<<pt_coord<<std::endl;
        }
    }
    //std::cout<<std::scientific<<std::setprecision(15)<<"H1 error for element "<<count<<" nu_err "<<nu_err<<" error "<<VMSerrH1<<std::endl;
    apf::setScalar(vmsErrH1,ent,0,VMSerrH1);

    VMSerrTotalL2 = VMSerrTotalL2+VMSerrL2*VMSerrL2;
    VMSerrTotalH1 = VMSerrTotalH1+VMSerrH1*VMSerrH1;

    apf::destroyElement(visc_elem);
    apf::destroyElement(pres_elem);
    apf::destroyElement(velo_elem);
    apf::destroyElement(velo_elem_old);
    apf::destroyElement(vof_elem);
    apf::destroyMeshElement(element);
    count++;
  } //end loop over elements
  
    PCU_Add_Doubles(&VMSerrTotalL2,1);
    PCU_Add_Doubles(&VMSerrTotalH1,1);
    if(PCU_Comm_Self()==0)
      std::cout<<std::scientific<<std::setprecision(15)<<"Total Error L2 "<<sqrt(VMSerrTotalL2)<<" H1 "<<sqrt(VMSerrTotalH1)<<std::endl;
    total_error = sqrt(VMSerrTotalH1);
    total_error_out = sqrt(VMSerrTotalH1);
    apf::destroyField(vmsErr);
    apf::destroyField(visc);

    getTotalLinearMomentum();
    apf::writeVtkFiles("CurrentErrorComputation",m);
    //std::abort();
} //end function

double get_nu_err(struct Inputs info){
  double stabTerm1 = 0.0;
  double stabTerm2 = 0.0;
  for(int i=0;i<info.nsd;i++){
    for(int j=0;j<info.nsd;j++){
       stabTerm1 += info.vel_vect[i]*info.gij[i][j]*info.vel_vect[j];
       stabTerm2 += info.gij[i][j]*info.gij[i][j];
     }
   } 
   double C_nu = 3.0;
   stabTerm2 = C_nu*info.visc_val*info.visc_val*sqrt(stabTerm2);
   double nu_err = 1.0/sqrt(info.visc_val*sqrt(stabTerm1) + stabTerm2);
   return nu_err;
}

apf::Vector3 getResidual(apf::Vector3 qpt,struct Inputs &info){
    apf::Vector3 grad_pres;
    apf::getGrad(info.pres_elem,qpt,grad_pres);
    //double visc_val = apf::getScalar(info.visc_elem,qpt);
    apf::Vector3 vel_vect;
    apf::getVector(info.velo_elem,qpt,vel_vect);
    apf::Vector3 vel_vect_old;
    apf::getVector(info.velo_elem_old,qpt,vel_vect_old);
    apf::Matrix3x3 grad_vel;
    apf::getVectorGrad(info.velo_elem,qpt,grad_vel);
    grad_vel = apf::transpose(grad_vel);

    apf::Vector3 tempConv;
    apf::Vector3 tempDiff;
    tempConv.zero();
    tempDiff.zero();
    for(int i=0;i<info.nsd;i++){
      for(int j=0;j<info.nsd;j++){
        tempConv[i] = tempConv[i] + vel_vect[j]*grad_vel[i][j];
      }
      tempConv[i] = tempConv[i] - info.g[i]; //body force contribution
    }
    apf::Vector3 tempResidual = (tempConv + grad_pres/info.density);
    //acceleration term
    tempResidual = tempResidual + (vel_vect-vel_vect_old)/dt_err;
    //std::cout<<"What is the acceleration contribution? "<<(vel_vect-vel_vect_old)/dt_err<<" vel_vect_old "<< vel_vect_old<<std::endl;

    info.vel_vect = vel_vect;
    info.grad_vel = grad_vel;
    //info.visc_val = visc_val;
    info.grad_pres = grad_pres;

  return tempResidual;
}

apf::Vector3 getResidual_debug(apf::Vector3 qpt,struct Inputs &info){
    apf::Vector3 grad_pres;
    apf::getGrad(info.pres_elem,qpt,grad_pres);
    //double visc_val = apf::getScalar(info.visc_elem,qpt);
    apf::Vector3 vel_vect;
    apf::getVector(info.velo_elem,qpt,vel_vect);
    apf::Vector3 vel_vect_old;
    apf::getVector(info.velo_elem_old,qpt,vel_vect_old);
    apf::Matrix3x3 grad_vel;
    apf::getVectorGrad(info.velo_elem,qpt,grad_vel);
    grad_vel = apf::transpose(grad_vel);

    apf::Vector3 tempConv;
    apf::Vector3 tempDiff;
    tempConv.zero();
    tempDiff.zero();
    apf::Vector3 convectiveTerm(0.0,0.0,0.0);
    for(int i=0;i<info.nsd;i++){
      for(int j=0;j<info.nsd;j++){
        tempConv[i] = tempConv[i] + vel_vect[j]*grad_vel[i][j];
        convectiveTerm[i] += vel_vect[j]*grad_vel[i][j];
      }
      tempConv[i] = tempConv[i] - info.g[i]; //body force contribution
    }
    apf::Vector3 tempResidual = (tempConv + grad_pres/info.density);
    //acceleration term
    tempResidual = tempResidual + (vel_vect-vel_vect_old)/dt_err;
    std::cout<<"What is the acceleration contribution? "<<(vel_vect-vel_vect_old)/dt_err<<" vel_vect_old "<< vel_vect_old<<std::endl;
    std::cout<<"grad pres "<<grad_pres/info.density<<std::endl;
    std::cout<<"gravity "<<info.g<<std::endl;
    std::cout<<"convective term "<<convectiveTerm<<std::endl;
    

    info.vel_vect = vel_vect;
    info.grad_vel = grad_vel;
    //info.visc_val = visc_val;
    info.grad_pres = grad_pres;

  return tempResidual;
}
