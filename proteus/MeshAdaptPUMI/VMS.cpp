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

inline void getProps(double*rho,double*nu)
//Function used to transfer MeshAdaptPUMIDrvr variables into global variables
{
  rho_0 = rho[0];
  nu_0 = nu[0];      
  rho_1 = rho[1];
  nu_1 = nu[1];
  return;
}

void MeshAdaptPUMIDrvr::get_VMS_error(double &total_error) 
{
  getProps(rho,nu);
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
  //*****               *****//

  //***** Compute the viscosity field *****//
  apf::Field* visc = getViscosityField(voff);

  apf::Field* vmsErr = err_reg = apf::createField(m,"VMSL2",apf::SCALAR,apf::getVoronoiShape(nsd,1));
  
  //Start computing element quantities
  int numqpt; //number of quadrature points
  int nshl; //number of local shape functions
  int elem_type; //what type of topology
  double weight; //value container for the weight at each qpt
  double Jdet;
  //apf::EntityShape* elem_shape;
  apf::Vector3 qpt; //container for quadrature points
  apf::MeshElement* element;
  apf::Element* visc_elem, *pres_elem,*velo_elem,*vof_elem;
  apf::Matrix3x3 J; //actual Jacobian matrix
  apf::Matrix3x3 invJ; //inverse of Jacobian
  apf::Matrix3x3 KJ(1.0,0.5,0.0,0.5,1.0,0.0,0.0,0.0,1.0); //Jacobian correction matrix for proper metric tensor
  apf::Matrix3x3 gij; //actual Jacobian matrix

  //apf::DynamicMatrix invJ_copy;
  //apf::NewArray <apf::DynamicVector> shdrv;
  //apf::NewArray <apf::DynamicVector> shgval_copy;
  
  apf::MeshIterator* iter = m->begin(nsd); //loop over elements
  apf::MeshEntity* ent;

  int count = 0 ; //for debugging purposes
  //Loop over elements and compute the VMS error in the L_2 norm
  while( (ent = m->iterate(iter)) ){

    element = apf::createMeshElement(m,ent);
    pres_elem = apf::createElement(pref,element);
    velo_elem = apf::createElement(velf,element);
    visc_elem = apf::createElement(visc,element);
    vof_elem = apf::createElement(voff,element);
    
    double strongResidualTauL2;
    double tau_m;
    double areaCheck=0.0;
    numqpt=apf::countIntPoints(element,int_order);

    
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
      
      if(count == 10){
        //std::cout<<"This is the tau "<<tau_m<<" velocity "<<vel_vect<<" grad vel "<<grad_vel<<" metric "<<gij<<" numdim "<<nsd<<" viscosity "<<visc_val<<" "<<KJ<<" gradP"<< grad_pres<<" stabTerms "<<stabTerm1<<" "<<stabTerm2<<std::endl;
        areaCheck += weight*Jdet;
      }

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
      if(count == 10 && k==0){
        std::cout<<" tempConv "<<tempConv<<" density "<<density<<" tempResidual "<<tempResidual<<" tempVal "<<tempVal<<std::endl;
      }

    } //end qpt loop

    double VMSerrL2 = sqrt(strongResidualTauL2);
    apf::setScalar(vmsErr,ent,0,VMSerrL2);
    if(count == 10){
      std::cout<<"This is the error "<< sqrt(strongResidualTauL2)<<std::endl;
      std::cout<<"This is the area check "<<areaCheck<<" true area is "<<apf::measure(m,ent);
    }
  
    apf::destroyElement(visc_elem);
    apf::destroyElement(pres_elem);
    apf::destroyElement(velo_elem);
    count++;
  } //end loop over elements

} //end function

  
