#include "MeshAdaptPUMI.h"
#include <PCU.h>
#include <petscksp.h>

#include <apf.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <apfDynamicMatrix.h>

#include <iostream>

enum {
  //PHI_IDX = 5
  VOF_IDX = 4
};

//used to attach error estimates to nodes
static void averageToEntity(apf::Field* ef, apf::Field* vf, apf::MeshEntity* ent) //taken from Dan's superconvergent patch recovery code
{
  apf::Mesh* m = apf::getMesh(ef);
  apf::Adjacent elements;
  m->getAdjacent(ent, m->getDimension(), elements);
  double s=0;
  for (std::size_t i=0; i < elements.getSize(); ++i)
    s += apf::getScalar(ef, elements[i], 0);
  s /= elements.getSize();
  apf::setScalar(vf, ent, 0, s);
  return;
}

static apf::Field* extractVOF(apf::Field* solution)
{
  apf::Mesh* m = apf::getMesh(solution);
  apf::Field* voff = apf::createLagrangeField(m,"proteus_vof",apf::SCALAR,1);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  apf::NewArray<double> tmp(apf::countComponents(solution));
  while ((v = m->iterate(it))) {
    apf::getComponents(solution, v, 0, &tmp[0]);
    double vof = tmp[VOF_IDX];
    apf::setScalar(voff, v, 0, vof);
  }
  m->end(it);
  return voff;
}

void getViscosity(double &nu_0, double &nu_1)
{
  //***** Viscosity hard-coded hook *****//
  nu_0 = 1.004e-6; //water 
  nu_1 = 1.5e-5; //air
  return;
}


void MeshAdaptPUMIDrvr::get_local_error() 
//This function aims to compute error at each element via ERM.
//First get the mesh and impose a 2nd order field
//Then get the desired quadrature points
{
  getViscosity(nu_0,nu_1);
  std::cout<<"Viscosity? "<<nu_0<<" "<<nu_1<<std::endl;

  //***** Get Phi First *****//
  apf::Field* voff = extractVOF(solution);
  apf::Field* visc = apf::createLagrangeField(m,"viscosity",apf::SCALAR,1);
  //*****               *****//
  
  //***** Compute the viscosity field *****//
  apf::Mesh*m = apf::getMesh(solution); 
  apf::MeshEntity* ent;
  apf::MeshIterator* iter = m->begin(0);
  double vof_val, visc_val;
  while(ent = m->iterate(iter)){ //loop through all elements
    vof_val=apf::getScalar(voff,ent,0);
    visc_val = nu_0*(1-vof_val)+nu_1*vof_val;
    apf::setScalar(visc, ent, 0,visc_val);
  }
  m->end(iter);
  //apf::writeVtkFiles("viscosity", m);

  //Start computing element quantities
  int numqpt; //number of quadrature points
  int nshl,nshl_sol; //number of local shape functions
  int elem_type; //what type of topology
  double weight; //value container for the weight at each qpt
  double Jdet;
  int approx_order = 2; //if using Lagrange shape functions. Serendipity is automatic
  apf::FieldShape* err_shape = apf::getLagrange(approx_order);
  apf::FieldShape* sol_shape = apf::getShape(solution);
  apf::EntityShape* elem_shape;
  apf::Vector3 qpt; //container for quadrature points
  apf::MeshElement* element;

  apf::Matrix3x3 J; //actual Jacobian matrix
  apf::Matrix3x3 invJ; //inverse of Jacobian
  apf::NewArray <double> shpval, shpval_sol; //array to store shape function values at quadrature points
  apf::NewArray <apf::Vector3> shgval, shgval_sol; //array to store shape function values at quadrature points

  apf::DynamicMatrix invJ_copy;
  apf::NewArray <apf::DynamicVector> shdrv,shdrv_sol;
  apf::NewArray <apf::DynamicVector> shgval_copy, shgval_sol_copy;
  
  iter = m->begin(m->getDimension()); //loop over elements
  //Vec F;
int testcount = 0;
  while(ent = m->iterate(iter)){ //loop through all elements
    element = apf::createMeshElement(m,ent);
    numqpt=apf::countIntPoints(element,approx_order);
    elem_type = m->getType(ent);
    nshl=apf::countElementNodes(err_shape,elem_type);
    nshl_sol=apf::countElementNodes(sol_shape,elem_type);
    if(testcount==0)
      std::cout<<"nshls "<<nshl<<" "<<nshl_sol<<std::endl;
    shpval.allocate(nshl);  shgval.allocate(nshl); shgval_copy.allocate(nshl); shdrv.allocate(nshl);
    shpval_sol.allocate(nshl_sol);  shgval_sol.allocate(nshl_sol); shgval_sol_copy.allocate(nshl_sol); shdrv_sol.allocate(nshl_sol);
    //loop through all qpts
    for(int k=0;k<numqpt;k++){
      apf::getIntPoint(element,approx_order,k,qpt); //get a quadrature point and store in qpt
      apf::getJacobian(element,qpt,J); //evaluate the Jacobian at the quadrature point
      if(J[2][2] == 0){ J[2][2]=1;} //can't compute an inverse with a 0 in a diagonal

      J = apf::transpose(J); //Is PUMI still defined in this way?
      invJ = invert(J);
      Jdet=apf::getJacobianDeterminant(J,m->getDimension()); 
      weight = apf::getIntWeight(element,approx_order,k);
      invJ_copy = apf::fromMatrix(invJ);

      //first get the shape function values for error shape functions
      elem_shape = err_shape->getEntityShape(elem_type);
      elem_shape->getValues(NULL,NULL,qpt,shpval);
      elem_shape->getLocalGradients(NULL,NULL,qpt,shgval); 

/*
      std::cout<<"Jdet "<<Jdet<<std::endl;
      std::cout<<"invJ "<<invJ_copy<<std::endl;
      std::cout<<"What is qpt and weight? "<<qpt<<" "<<weight<<std::endl;
*/
      for(int i =0;i<nshl;i++){ //get the true derivative
        shgval_copy[i] = apf::fromVector(shgval[i]);
        //std::cout<<"What is shgval? "<<shgval_copy[i]<<std::endl;
        apf::multiply(shgval_copy[i],invJ_copy,shdrv[i]); 
//        if(testcount==0 && k==0)
//        std::cout<<"What is shdrv? "<<shdrv[i]<<std::endl;
      }

      //solution shape functions
      elem_shape = sol_shape->getEntityShape(elem_type);
      elem_shape->getValues(NULL,NULL,qpt,shpval_sol);
      elem_shape->getLocalGradients(NULL,NULL,qpt,shgval_sol); 

      for(int i =0;i<nshl_sol;i++){ //get the true derivative
        shgval_sol_copy[i] = apf::fromVector(shgval_sol[i]);
        //std::cout<<"What is shgval? "<<shgval_copy[i]<<std::endl;
        apf::multiply(shgval_sol_copy[i],invJ_copy,shdrv_sol[i]); 
//        if(testcount==0 && k==0)
//        std::cout<<"What is shdrv_sol? "<<shdrv_sol[i]<<std::endl;
      }


  //Store into matrix

/*
    VecCreate(PETSC_COMM_SELF,&F); 
    VecSetSizes(F,5,5);
    VecSetFromOptions(F); 
    VecSetUp(F);
    for(int i=0;i<5;i++){
      VecSetValue(F,i,2*i+1,INSERT_VALUES);
    }
    VecView(F,PETSC_VIEWER_STDOUT_SELF);
    VecDestroy(&F);

*/
    }

testcount++;
  }
  m->end(iter);
  apf::destroyField(voff);
  apf::destroyField(visc);
  printf("It cleared the function.\n");
}


